library(tidyverse)
library(raster)
library(sf)
library(broom)
library(RColorBrewer)
library(classInt)
library(gridExtra)
library(tmap)

calname = "Sacks_ZARC_fill_fill_120d"
# ddversionstring = "bcagcfsr_test"
# ddversionstring = "agcfsr_test_small"
ddversionstring = "rgagcfsr"

betasfolder = paste0("gdd_betas_pixel_partial/",calname,"/",ddversionstring,"/")
# infolder  = "xavier_computed/"
# outfolder  = "xavier_dfs/"

# dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# crops = c("Soybeans")
# years = 2002:2008

crops = c("Maize","Soybeans","Cotton")
years = 2002:2008
# years = 2002:2003
# deltats = 0:5 # Must get zero

# Each crop will use different tresholds for GDD and EDD.
# Give lower limit EDD names for GDD, and later GDD can be calculated as EDDl-EDDu
# Values taken from Schlenker and Roberts can be found in digitized/traced.ods
extvnames = c("edd38","edd38","edd38")
eddvnames = c("edd29","edd30","edd32")
gddvnames = c("edd10","edd10","edd15")
names(extvnames) <- crops
names(eddvnames) <- crops
names(gddvnames) <- crops

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

# Points shapefile
shpfname = 'GIS/COLROW30.shp'
cntfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'

# Reference raster
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

# SR values to compare to
srrefdata = read.csv("aux_figures/test_table_A5.csv") %>% filter(source == "SR Table A5")

# Function definitions
# Evaluates DD betas in a data.frame for a given deltat
eval_dd <- function(data,ddname,deltat) {
  out = data[paste0(ddname,"beta1")]*deltat +
    data[paste0(ddname,"beta2")]*(deltat^2) +
    data[paste0(ddname,"beta3")]*(deltat^3) 
  names(out) <- paste0("d",ddname)
  return(out)
}
# Evaluates PARTIAL method dyld (fractional) in a data frame given sensitivities and a fixed scale factor
# Must decide whether to use fgff or cgdd
eval_partial_dyld_fixscale <- function(data,cropeddsens,cropgddsens,scalesens,deltat) {
  dyld = eval_dd(data, "agdd", deltat)*cropgddsens + 
    eval_dd(data, "bgdd", deltat)*cropeddsens +
    eval_dd(data, "fgdd", deltat)*cropeddsens*scalesens
  
  # dyld = eval_dd(data, "gdd", deltat)*cropgddsens +
  # eval_dd(data, "edd", deltat)*cropeddsens*scalesens
  names(dyld) <- "dyld"
  return(dyld)
}
# Evaluates dyld (fractional) in a data frame given sensitivities and a fixed scale factor
eval_dyld_fixscale <- function(data,cropeddsens,cropgddsens,scalesens,deltat) {
  dyld = eval_dd(data, "gdd", deltat)*cropgddsens +
    eval_dd(data, "edd", deltat)*cropeddsens*scalesens
  names(dyld) <- "dyld"
  return(dyld)
}
# Evaluates dyld (fractional) in a data frame given sensitivities and a stepwise scale factor
eval_dyld_stepscale <- function(data,cropeddsens,cropgddsens,deltat,minscale,slp,inter) {
  scalesens = inter + slp*data$tmaxbase
  scalesens = pmin(1,scalesens)
  scalesens = pmax(minscale[1],scalesens)
  dyld = eval_dd(data, "gdd", deltat)*cropgddsens +
    eval_dd(data, "edd", deltat)*cropeddsens*scalesens
  names(dyld) <- "dyld"
  return(dyld)
}


eval_dyld_countries_fixscale <- function(data,countries,cropeddsens,cropgddsens,scalesens,deltat) {
  filtdata = data %>% filter(ADM0_A3 %in% countries)
  filtdata$dyld = eval_dyld_fixscale(filtdata,cropeddsens,cropgddsens,scalesens,deltat)
  filtdata %>% na.omit() %>%
    group_by(ADM0_A3) %>% 
    summarise(dyld = sum(dyld*area)/sum(area)) %>% 
    return()
}

eval_dyld_countries_stepscale <- function(data,countries,cropeddsens,cropgddsens,deltat,minscale,slp,inter) {
  filtdata = data %>% filter(ADM0_A3 %in% countries)
  filtdata$dyld = eval_dyld_stepscale(filtdata,cropeddsens,cropgddsens,deltat,minscale,slp,inter)
  filtdata %>% na.omit() %>%
    group_by(ADM0_A3) %>% 
    summarise(dyld = sum(dyld*area)/sum(area)) %>% 
    return()
}

# deltat = 5
# data = cropdata
# countries = c("USA","BRA")
# eval_dyld_countries_fixscale(cropdata,countries,eddsens[[crop]],gddsens[[crop]],1.0,deltat)

# Reading shapefiles and creating country table
shp = st_read(shpfname)
cntshp = st_read(cntfname)
st_crs(shp) <- st_crs(cntshp)
cntshp = st_join(shp,cntshp)[c("COLROW30","ADM0_A3")]
cntdata = cntshp %>% st_set_geometry(NULL)

# Reference raster
refrast = raster(reffname)

# Base temperatures at the GLOBIOM grid
tempbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",calname,"/Tbaseline_tmp.csv")
tmaxbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",calname,"/Tbaseline_tmx.csv")

tempbase = read.csv(tempbasefname) %>% gather("Crop","tempbase",-ID)
tmaxbase = read.csv(tmaxbasefname) %>% gather("Crop","tmaxbase",-ID)
tempdata = left_join(tmaxbase,tempbase, by = c("ID","Crop")) %>%
  dplyr::rename(COLROW30 = "ID")

areadata = read.csv("GIS/areas.csv")

repdata = data.frame()
brrepdata = data.frame()

# crop = "Soybeans"
for (crop in crops) {
# for (crop in c("Soybeans","Maize")) {
  # Read betas for that crop
  vnames = c("agdd","bgdd","cgdd","fgdd")
  for (vname in vnames) {
    betas = read.csv(paste0(betasfolder,"/",crop,".betas.",vname,".csv"))
    names(betas)[-1] <- paste0(vname,names(betas)[-1])
    if (vname == vnames[1]) {
      allbetas = betas
    } else {
      allbetas = left_join(allbetas,betas)
    }
  }
  
  
  cropareadata = areadata[c("COLROW30",crop)] %>% dplyr::rename(area = crop)
  cropdata = cntdata %>% 
    left_join(cropareadata) %>% 
    left_join(allbetas) %>% 
    left_join(filter(tempdata,Crop == crop)) %>%
    mutate(COLROW30 = as.factor(COLROW30))
  
  # data = filter(cropdata,COLROW30 %in% c("252 - 212", "253 - 212"))
  usdata = filter(cropdata, ADM0_A3 == "USA")
  brdata = filter(cropdata, ADM0_A3 == "BRA")
  
  
  for (deltat in 0:10) {
    usdata[paste0("dyld",deltat)] = eval_partial_dyld_fixscale(usdata,eddsens[[crop]],gddsens[[crop]],1.0,deltat)
    
    out = usdata %>% na.omit %>% summarise(poi = weighted.mean(.[[paste0("dyld",deltat)]],.$area)) %>% as.numeric()
    repdata = rbind(repdata, data.frame(deltat = deltat, Crop = crop, dyld = out*100, source = "New Approach"))
    
    brdata[paste0("dyld",deltat)] = eval_partial_dyld_fixscale(brdata,eddsens[[crop]],gddsens[[crop]],1.0,deltat)
    
    out = brdata %>% na.omit %>% summarise(poi = weighted.mean(.[[paste0("dyld",deltat)]],.$area)) %>% as.numeric()
    brrepdata = rbind(brrepdata, data.frame(deltat = deltat, Crop = crop, dyld = out*100, source = "New Approach"))
  }
  
  if (iplotglobal) {
    deltats = seq(1,6)
    scennames = paste0("deltaT",sprintf("%1i",deltats))
    for (i in 1:length(deltats)) {
      cropdata[scennames[i]] = eval_partial_dyld_fixscale(cropdata,eddsens[[crop]],gddsens[[crop]],1.0,deltats[i])*100.0
    }
    cropshp = left_join(shp,filter(cropdata,Crop == crop),by = c("COLROW30" = "COLROW30"))
    valid = scennames
    croprast = rasterize(cropshp, refrast, field = valid, fun = mean, background = NA_real_,
                         by = NULL)
    # croprasts[[crop]] <- croprast
    
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
  }
  
}
allrepdata = rbind(srrefdata,repdata)

ggplot(allrepdata) + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source),size=2) +
  ggtitle("Comparison with U.S. average impacts", subtitle = ddversionstring) + 
  xlim(0,6) +
  ylim(-100,0)

ggplot(brrepdata) + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source),size=2) +
  ggtitle("Brazil average impacts", subtitle = ddversionstring) + 
  xlim(0,6) +
  ylim(-100,0)

bothdata = rbind(mutate(allrepdata, Region = "USA"),
  mutate(brrepdata, Region = "Brazil"))
ggplot(filter(bothdata, source == "New Approach")) + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = Region),size=2) +
  xlim(0,6) +
  ylim(-120,0)



ggplot(filter(allrepdata, Crop == "Soybeans")) + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source),size=2) +
  xlim(0,6) +
  ylim(-150,0)

usshpdata = cntshp %>% filter(ADM0_A3 == "USA") %>% left_join(usdata)

plot(usshpdata[c("dyld5","geometry")], xlim = c(-100,-90), pch = 15)


plot(usshpdata[c("dyld6","geometry")])


eval_dyld_countries_stepscale(data,c("USA","BRA"),eddsens[[crop]],gddsens[[crop]],1.0,5,minscale,slp,inter)


eval_dyld_stepscale(data,eddsens[[crop]],gddsens[[crop]],1.0,1,minscale,slp,inter)

poi = data.frame()
for (deltat in 0:6) {
  # dum = eval_dyld_countries_fixscale(data,c("USA","BRA"),eddsens[[crop]],gddsens[[crop]],1.0,deltat)
  dum = eval_dyld_countries_stepscale(cropdata,c("USA","BRA"),eddsens[[crop]],gddsens[[crop]],deltat,minscale,slp,inter)
  dum$deltat = deltat
  poi = rbind(poi,dum)
}

ggplot(poi) + geom_line(aes(x=deltat,y=dyld,color=ADM0_A3), size = 2) +
  ggtitle("Same scaling as previous approach (Soybens)")

sminscale = 0.75
sslp = -0.1
sinter = 3.8

f1(sminscale,-0.4,sinter)


f1 = function(minscale,slp,inter) {
  dum = eval_dyld_countries_stepscale(cropdata,c("USA","BRA"),eddsens[[crop]],gddsens[[crop]],3,minscale,slp,inter)
  usval = as.numeric(dum[dum$ADM0_A3 == "USA","dyld"])
  brval = as.numeric(dum[dum$ADM0_A3 == "BRA","dyld"])
  return(c(usval,brval,brval/usval))
}

pltfun = function(minscale,slp,inter) {
  pdata = data.frame(tmax = seq(15,40,0.1))
  pdata$scalesens = inter + slp*pdata$tmax
  pdata$scalesens = pmin(1,pdata$scalesens)
  pdata$scalesens = pmax(0,pdata$scalesens)
  ggplot(pdata) + geom_line(aes(x=tmax,y=scalesens)) + ylim(0,1)
}

pltfun(sminscale,-0.3,5)

lotmax = 20
hitmax = 32
minscale = 0.5
fit=lm(scalesens~tmax,data.frame(tmax=c(lotmax,hitmax),scalesens=c(1,minscale)))
f1(minscale,fit$coefficients[2],fit$coefficients[1])



# deltats = seq(1,6)
# scennames = paste0("deltaT",sprintf("%1i",deltats))
# for (i in 1:length(deltats)) {
#   cropdata[scennames[i]] = eval_partial_dyld_fixscale(cropdata,eddsens[[crop]],gddsens[[crop]],1.0,deltats[i])*100.0
# }
# cropshp = left_join(shp,filter(cropdata,Crop == crop),by = c("COLROW30" = "COLROW30"))
# valid = scennames
# croprast = rasterize(cropshp, refrast, field = valid, fun = mean, background = NA_real_,
#                      by = NULL)
# # croprasts[[crop]] <- croprast
# 
# breaks = seq(-100,100,20)
# pal = brewer.pal(n = length(breaks), name = "RdBu")
# tit = crop
# 
# plotobj <- tm_shape(croprast, bbox = c(-130,160,-50,65)) + 
#   tm_raster(palette = pal, breaks = breaks,
#             title = expression(paste(Delta,"Yield (%)"))) + 
#   tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
#             legend.outside = T) + 
#   tm_layout(panel.show = T, panel.labels = scennames,
#             panel.label.size = 1.2,
#             main.title.position = "center", main.title = tit, main.title.size = 1.0) +
#   tm_facets(nrow = 3, ncol = 2)
# print(plotobj)
