Packages <- c("Rcpp","dplyr","tidyverse","gridExtra","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# Function that applies the DD~T curves
source("eval_nxdd.R")

tdbase = "Sacks_ZARC_fill_fill_120d"

ddversionstring = "agcfsr_unbound_usa"
# ddversionstring = "agcfsr_bound_usa"
# ddversionstring = "merra"
# ddversionstring = "wgt_crop"
# ddversionstring = "wgt_crop_un1bound"

crops = c("Maize","Soybeans","Cotton")
# crops = c("Maize","Soybeans")

# Base sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Reference countries and what impacts should be at hirefcountry's Tmax at reference deltat. 
hirefcountry = "BRA"
lorefcountry = "USA"
refdeltat = 6.0
# refmult = 1.25
refmult = 1.50
refvals = c(-55.88, -48.82, -46.72)*refmult
names(refvals) <- crops
# Use computed instead of the reported refvals
icomputerefvals = TRUE

# XDD ~ Tx functions
betasfolder = paste0("gdd_betas/",tdbase,"/",ddversionstring,"/")

# Base data to evaluate functions on (tmx and tmp are Tbaselines)
bigbasedata = readRDS("tbaseline_countries.rds")

# Output file name
# outfolder = paste0("tmaxsens/",ddversionstring,"/")
# outfname = paste0(outfolder,"tmaxsens_linear_2_",refmult,".csv")

# SR values to compare to
srrefdata = read.csv("aux_figures/test_table_A5.csv") %>% filter(source == "SR Table A5")

# GLOBIOM grid, points shapefile and reference for rasterization
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
countriesfname = "GIS/ne_50m_admin_0_countries_lakes.shp"
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"
areafname = "GIS/areas.csv"

shp = st_read(shpfname)
refrast = raster(reffname)
convdata = read.csv(convfname) %>% rename(GLOBIOM_ID = "final.ID")
shp <- left_join(shp,convdata, by = c("COLROW30" = "final.COLROW30"))

areadata = read.csv(areafname) %>% gather(Crop,areacrop,-ID,-Country)

# Also read countries and add both codes and English names as variables 
countriesshp = st_read(countriesfname)
countriesshp %>% dplyr::select(ADM0_A3,NAME_EN,geometry) ->countriesshp
st_crs(shp) <- st_crs(countriesshp) #FIXME: This assumes shp's SRS is undefined, but is WGS84
st_join(shp,countriesshp) -> shp

# Plot maps of impacts on the GLOBIOM grid
iglobplot = TRUE

# Base temperatures at the GLOBIOM grid
tempbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmp.csv")
tmaxbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmx.csv")

tempbase = read.csv(tempbasefname) %>% gather("Crop","tempbase",-ID)
tmaxbase = read.csv(tmaxbasefname) %>% gather("Crop","tmaxbase",-ID)
tempdata = left_join(tmaxbase,tempbase, by = c("ID","Crop"))


# dir.create(outfolder,showWarnings = F, recursive = T)

# caldata = data.frame(tmx = c(27,27,30), X.Temp = c(1,6,1), deltay = c(-6.4,-55.9,-1))
# caldata = read.csv("../moore_coefficients/includedstudies_revisions_tbaseline.csv")
# names(caldata)[names(caldata) == "Yield.change...."] <- "deltay"


# Functions to evaluate the method. Depend on a lot of variables being in .GlobalEnv
# Evaluate using a fixed scalesens. Set it to 1 to get the original value
eval_fixscale_point <- function(tmx,tmp,deltat,scalesens) {
  # tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
  # return(tmp)
}
# Evaluate using a stepwise tmaxbase~scalesens relationship
eval_scaled_step_point <- function(tmx,deltat,coeftab) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  scalesens = coeftab$inter[1] + coeftab$slp[1]*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(coeftab$minscale[1],scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}

# Evaluate using a stepwise tmaxbase~scalesens relationship, taking tmp as well
eval_scaled_step_point_tmp <- function(tmx,tmp,deltat,cropcoeftab) {
  # tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  scalesens = cropcoeftab$inter[1] + cropcoeftab$slp[1]*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(cropcoeftab$minscale[1],scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}


coeftab = data.frame()
plots_range = list()
plots_fun = list()
croprasts = list()
simdata = data.frame()
# crop = "Maize"
for (crop in crops) {
  
  # Read the crop-specific coefficients for the T-GDD functions
  betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
  betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))
  
  cropshp <- shp %>%
    left_join(filter(tempdata,Crop == crop),by = c("COLROW30" = "ID")) %>%
    left_join(filter(areadata,Crop == crop),by = c("GLOBIOM_ID" = "ID","Crop" = "Crop")) %>%
    mutate(areacrop = ifelse(is.na(areacrop),0,areacrop))
  deltats = seq(1,6)
  scennames = paste0("deltaT",sprintf("%1i",deltats))
  scentable = data.frame(scenname = scennames, deltat = deltats)
  for (i in 1:length(deltats)) {
    # cropshp[scennames[i]] = eval_scaled_step_point_tmp(cropshp$tmaxbase,cropshp$tempbase,deltats[i],cropcoeftab)*100.0
    cropshp[scennames[i]] = eval_fixscale_point(cropshp$tmaxbase,cropshp$tempbase,deltats[i],1.0)*100.0
  }
  
  cropdata = cropshp %>% gather(deltat,dyld,scennames) %>% 
    mutate(deltat = as.integer(substr(deltat,7,10))) %>%
    st_set_geometry(NULL)
  
  cropdata %>% filter(ADM0_A3 == "USA") %>% 
    group_by(ADM0_A3,deltat,Crop) %>% 
    summarise(aveimpact = weighted.mean(dyld,areacrop,na.rm = T)) %>%
    ungroup %>% mutate(source = "Impact average(USA)") %>% rbind(simdata,.) ->simdata
  
  cropdata %>% filter(ADM0_A3 == "USA") %>% filter(XCOORD >= -100.0) %>% 
    group_by(ADM0_A3,deltat,Crop) %>% 
    summarise(aveimpact = weighted.mean(dyld,areacrop,na.rm = T)) %>%
    ungroup %>% mutate(source = "Impact average(USA East)") %>% rbind(simdata,.) ->simdata
  
  basedata = filter(bigbasedata,Crop == crop)
  ustmx = filter(basedata, country == "USA")[["tmx"]]
  ustmp = filter(basedata, country == "USA")[["tmx"]]
  avsimdata = data.frame(ADM0_A3 = "USA", deltat = deltats, Crop = crop,
                        aveimpact = eval_fixscale_point(ustmx,ustmp,deltats,1.0)*100.0,
                        source = "Average Temp (USA)")
  simdata = rbind(simdata,avsimdata)
  
  
  valid = scennames
  croprast = rasterize(cropshp, refrast, field = valid, fun = mean, background = NA_real_,
                       by = NULL)
  croprasts[[crop]] <- croprast
  
  breaks = seq(-100,100,20)
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  tit = crop
  
  plotobj <- tm_shape(croprast, bbox = c(-130,160,-50,65)) +
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5,
              legend.outside = T) +
    tm_layout(panel.show = T, panel.labels = scennames,
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = 3, ncol = 2)
  print(plotobj)
  
  
  # Synthetic data for plotting 
  pdata = data.frame(tmx = seq(15,36,3))
  pdata$tmp = fittmx$coefficients[[1]] + pdata$tmx*fittmx$coefficients[[2]]
  dumdata = pdata
  pdata = data.frame()
  for (dt in seq(0,6,0.5)) {
    dumdata$deltat = dt
    pdata = rbind(pdata,dumdata)
  }
  
  # # Plot the estimated range of impacts
  # plotrange = ggplot(pdata) + 
  #   geom_line(aes(x = deltat, y = dyld*100, color = as.factor(tmx)), size = 1.5) + 
  #   xlab(expression(paste(Delta,"T (°C)"))) +
  #   ylab(expression(paste(Delta,"T (pp)"))) +
  #   scale_colour_discrete(name="Tmax") +
  #   labs(title = crop) + ylim(-80,50)
  
  
  # plots_range[[length(plots_range)+1]] = plotrange
  

  
}

simdata %>% select(-ADM0_A3) %>% rename(dyld = aveimpact) %>%
  rbind(srrefdata,.) -> allsimdata



# Plot of composite C3 crops
if (iglobplot) {
  breaks = seq(-100,100,20)
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  tit = "Composite C3 (average of Soybeans and Cotton)"
  
  croprast = (croprasts$Soybeans + croprasts$Cotton )/2.0
  
  plotobj <- tm_shape(croprast, bbox = c(-130,160,-50,65)) + 
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = scennames,
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = 3, ncol = 2)
  print(plotobj)
  
  write.csv(coeftab, file = outfname, row.names = F)
}

ggplot(allsimdata) + 
  geom_line(aes(x = deltat, y = dyld, color = source), size = 2) +
  facet_wrap(~Crop) +
  xlim(1,6) + ylim(-60,10)

# Zoom in USA Cotton
breaks = seq(-40,40,10)
pal = brewer.pal(n = length(breaks), name = "RdBu")
tit = crop
usreg = st_read("GIS/states_us_conus.shp")

plotobj <- tm_shape(croprasts[[crop]], bbox = usreg) +
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Yield (%)"))) +
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5,
            legend.outside = T) +
  tm_shape(usreg) + tm_borders() +
  tm_layout(panel.show = T, panel.labels = scennames,
            panel.label.size = 1.2,
            main.title.position = "center", main.title = tit, main.title.size = 1.0) +
  tm_facets(nrow = 3, ncol = 2)
print(plotobj)



# plots_all = list()
# for (i in seq(1,(2*length(plots_range)),2)) {
#   plots_all[i] = plot
# }

do.call(grid.arrange,c(plots_range,plots_fun,c(nrow = 2)))

ggplot(srrefdata) +
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source), size = 1.5)

# grid.arrange(plotrange,plotfun,nrow = 1)
