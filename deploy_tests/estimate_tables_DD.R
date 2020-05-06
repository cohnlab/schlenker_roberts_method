Packages <- c("Rcpp","ggridges","gridExtra","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

infolder = "deploy_tests/Coefficient Tables/"
outfname = "deploy_tests/results_table_DD.rds"

# Tbaselines
tempbasefname = paste0(infolder,"DD approach/Tbaseline_tmp.csv")
tmaxbasefname = paste0(infolder,"DD approach/Tbaseline_tmx.csv")

# Agroserv-Temp coefficients
nlcoeffname = paste0(infolder,"AgroservTemp/NL.csv")
lcoeffname = paste0(infolder,"AgroservTemp/L.csv")

# Folder with the DD~T tables
betasfolder = paste0(infolder,"DD approach/")

# EDD_sensitivity_scaling_factor
scalefacfname = paste0(infolder,"DD approach/EDD_sensitivity_scaling_factor.csv")

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

# Wrapper for cpp_eval_nxdd. Forces a single-column data.frame into a vector
eval_nxdd <- function(betas,tvec) {
  if (is.data.frame(tvec)) {
    tvec = tvec[,]
  }
  # Using Rcpp here to avoid figuring out a faster lookup in R
  # Evaluete a nEDD or nGDD function at a given T(max,mean)
  cppFunction('NumericVector cpp_eval_nxdd(NumericMatrix x, NumericVector tvec) {
  int nrow = x.nrow(), ncol = x.ncol();
  int vsize = tvec.size();
  int p = 0;
  double t;
  double val;
  NumericVector vout(vsize);
  for (int p = 0; p < vsize; ++p) {
  t = tvec(p);
  for (int i = 0; i < nrow; ++i) {
  if (t >= x(i,0) && t<= x(i,1)) {
    vout(p) = x(i,2) + x(i,3)*t;
    break;
  }
  }
  }
  return vout;
}')
  return(cpp_eval_nxdd(as.matrix(betas),as.vector(tvec)))
}

# Evaluate DD method using a stepwise tmaxbase~scalesens relationship, taking tmp as well
eval_dd_scaled_step_point_tmp <- function(tmx,tmp,deltatmx,deltatmp,cropcoeftab) {
  # tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltatmp) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltatmx) - eval_nxdd(betasedd,tmx)
  
  scalesens = cropcoeftab$inter[1] + cropcoeftab$slp[1]*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(cropcoeftab$minscale[1],scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
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
    betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
    betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))
    
    cropcoeftab = read.csv(scalefacfname) %>% filter(Crop == crop)
    
    
    cropdata$dyld = with(cropdata,
                         eval_dd_scaled_step_point_tmp(tmaxbase,tempbase,deltatmax,deltatemp,cropcoeftab) )
    
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
