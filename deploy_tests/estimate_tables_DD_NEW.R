Packages <- c("Rcpp","ggridges","gridExtra","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

infolder = "deploy_tests/Coefficient Tables/"
outfname = "deploy_tests/results_table_DD_NEW.rds"

# Tbaselines
tempbasefname = paste0(infolder,"DD approach/Tbaseline_tmp.csv")
tmaxbasefname = paste0(infolder,"DD approach/Tbaseline_tmx.csv")

# Agroserv-Temp coefficients
nlcoeffname = paste0(infolder,"AgroservTemp/NL.csv")
lcoeffname = paste0(infolder,"AgroservTemp/L.csv")

# Folder with the DD~T tables
# betasfolder = paste0(infolder,"DD approach/")
betasfolder = paste0("gdd_betas_pixel_fill/Sacks_ZARC_fill_fill_120d/rgagcfsr_bound_1995_2005/")

# EDD_sensitivity_scaling_factor
# scalefacfname = paste0(infolder,"DD approach/EDD_sensitivity_scaling_factor.csv")

# Base sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
crops = c("Maize","Soybeans")
eddsens = c(-0.006435774107513,-0.005894686874659)
gddsens = c(0.000317954802039,0.000403658172117)
names(eddsens) <- crops
names(gddsens) <- crops

# This is the file with GHG deltaTs for the only scenario we're evaluating
# Tables for the other scenarios can be built on this one
# We'll also use only one year
# rcpdeltatfname = paste0(infolder,"AgroservTemp/deltaTssp245.csv")
rcpdeltatfname = paste0(infolder,"AgroservTemp/deltaTrcp245.csv")
rcpscenyear = 2050
rcpscenname = "ssp245"

dluLs = c(0.00,0.25,0.50,0.75,1.00)
dluNLs = c(0.00,0.25,0.50,0.75,1.00)

scens = c(rcpscenname,"deltaT0","deltaT1","deltaT2","deltaT3","deltaT4","deltaT5")

###### Function definitions

# Evaluates DD betas in a data.frame for a given deltat
eval_dd <- 
  function(data,ddname,deltat) {
    out = data[paste0(ddname,"beta1")]*deltat +
      data[paste0(ddname,"beta2")]*(deltat^2) +
      data[paste0(ddname,"beta3")]*(deltat^3) 
    names(out) <- paste0("d",ddname)
    return(out)
  }
# Evaluates dyld (fractional) in a data frame given sensitivities and a fixed scale factor
# Must decide whether to use fgff or cgdd
eval_dyld_new <- 
  function(data,cropeddsens,cropgddsens,deltat) {

    dyld = eval_dd(data, "GDD", deltat)*cropgddsens +
      eval_dd(data, "EDD", deltat)*cropeddsens
    dyld = exp(dyld)-1
    dyld = as.numeric(dyld[,])
    return(dyld)
  }

# Begin script

# First we'll make a table with all input data that is not scenario dependent
tempbase = read.csv(tempbasefname) %>% gather("Crop","tempbase",-ID) #FIXME there one bloody NA here
tmaxbase = read.csv(tmaxbasefname) %>% gather("Crop","tmaxbase",-ID)
tempdata = left_join(tmaxbase,tempbase, by = c("ID","Crop"))

nlcoefdata = read.csv(nlcoeffname)
lcoefdata = read.csv(lcoeffname) 
agtdata = rbind(nlcoefdata,lcoefdata) %>% spread("PARAMETER","Value")

basedata = left_join(tempdata, agtdata, by = c("ID"))
basedata$Crop = as.factor(basedata$Crop)

# Now build a version with all dlu scenarios
# FIXME: Looping here is not very smart, some weird join might be faster
dlubasedata = data.frame()
for (dluL in dluLs) {
  for (dluNL in dluNLs) {
    tbasedata = basedata
    tbasedata$dluL = dluL
    tbasedata$dluNL = dluNL
    dlubasedata = rbind(dlubasedata,tbasedata)
  }
}

# Now get agtdeltat, the warming due to dlus
dlubasedata$agtdeltatmax = with(dlubasedata, NL*dluNL + L*dluL)
dlubasedata$agtdeltatemp = with(dlubasedata, NL*dluNL/2.0 + L*dluL/2.0)

# Read the RCP deltat table, but we won't do anything with it until looping scenarios
rcpdeltatdata = read.csv(rcpdeltatfname) %>% 
  filter(ScenYear == rcpscenyear) %>% 
  spread(PARAMETER,Value) %>% 
  dplyr::select(-scen,-ScenYear) # Remove the scen variable, we'll build our own

# BEGIN SCENARIO LOOP
outalldata = data.frame()
for (scen in scens) {
  print(scen)
  # Here we join if the scenario is the RCP one, or get the fixed deltaT from the last 
  # character of scenname 
  if (scen == rcpscenname) {
    scendata = dlubasedata %>% left_join(rcpdeltatdata, by = c("ID","Crop"))  
  } else {
    deltat = as.numeric(substr(scen,nchar(scen),nchar(scen)))
    scendata = dlubasedata
    scendata$gccdtemp = deltat
    scendata$gccdtmax = deltat
  }
  scendata$scen = scen
  
  # Now we evalueate actual deltats
  scendata$deltatemp = scendata$agtdeltatemp + scendata$gccdtemp
  scendata$deltatmax = scendata$agtdeltatmax + scendata$gccdtmax
  
  outscendata = data.frame()
  # BEGIN CROP LOOP
  for (crop in crops) {
    cropdata = scendata %>% filter(Crop == crop)
    
    # Read the crop-specific coefficients for the T-GDD functions
    vnames = c("GDD","EDD")
    for (vname in vnames) {
      # betas = read.csv(paste0(betasfolder,"/",crop,".betas.",vname,".csv"))
      betas = read.csv(paste0(betasfolder,"/",crop,".betas",vname,".csv"))
      names(betas)[-1] <- paste0(vname,names(betas)[-1])
      if (vname == vnames[1]) {
        allbetas = betas
      } else {
        allbetas = left_join(allbetas,betas)
      }
    }
    allbetas = allbetas %>% rename (ID = "COLROW30")
    cropdata = cropdata %>% left_join(allbetas)
    
    # cropcoeftab = read.csv(scalefacfname) %>% filter(Crop == crop)
    
    cropdata$dyld = eval_dyld_new(cropdata,eddsens[[crop]],gddsens[[crop]],cropdata$deltatemp) %>%
      as.numeric()
    
    outscendata =
      cropdata %>% 
      dplyr::rename(dtemp = deltatemp, dtmax = deltatmax) %>%
      dplyr::select(ID,scen,Crop,dluL,dluNL,dyld,dtemp,dtmax) %>%
      bind_rows(outscendata)
  }
  outalldata = outscendata %>% bind_rows(outalldata)
}

# Convert factors
outalldata$scen <- as.factor(outalldata$scen)
outalldata$Crop <- as.factor(outalldata$Crop)
outalldata$dluL <- as.factor(outalldata$dluL)
outalldata$dluNL <- as.factor(outalldata$dluNL)

# dyld into percentages
outalldata$dyld = outalldata$dyld*100.0

saveRDS(outalldata,outfname)
