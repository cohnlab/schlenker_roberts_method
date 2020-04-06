Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lspline)
library(classInt)

#Here we estimate the DD~T functions and apply SR to US average values

calname = "Sacks_ZARC_fill_fill_120d"
# ddversionstring = "merra"
# ddversionstring = "merra_bound"
# ddversionstring = "agcfsr_bound"
# ddversionstring = "agcfsr_bound_usa"
ddversionstring = "agcfsr_unbound_usa"
# ddversionstring = "agmerra_bound"
# ddversionstring = "agmerra_bound_usa"
# ddversionstring = "wgt_crop_unbound"

# infolder  = paste0("sheffield_dfs_vpd/",calname,"/","wgt_crop","/")
# infolder  = paste0("merra_dfs/",calname,"/","merra","/")
infolder  = paste0("agcfsr_dfs/",calname,"/","agcfsr","/")
# infolder  = paste0("agmerra_dfs/",calname,"/","agmerra","/")
outfolder = paste0("gdd_betas/",calname,"/",ddversionstring,"/")

# years = 2009:2014
years = 2002:2008
# years = 1991:2008
# years = c(2002)
# years = c(2003)
# years = 2002:2003

# Filter specific countries?
ifiltcnt = TRUE
filtcnts = c("USA")

# crops = c("Soybeans")
crops = c("Maize","Soybeans","Cotton")

# Each crop will use different tresholds for GDD and EDD.
# Give lower limit EDD names for GDD, and later GDD can be calculated as EDDl-EDDu
# Values taken from Schlenker and Roberts can be found in digitized/traced.ods
extvnames = c("edd38","edd38","edd38")
eddvnames = c("edd29","edd30","edd32")
gddvnames = c("edd10","edd10","edd15")
names(extvnames) <- crops
names(eddvnames) <- crops
names(gddvnames) <- crops
# Should EDDs have an upper bound?
ieddbound = FALSE

# Cuts (hinges, breaks) of the spline approximation
cuts = seq(9,37,2)

####### SR stuff

# SR values to compare to
srrefdata = read.csv("aux_figures/test_table_A5.csv") %>% filter(source == "SR Table A5")

# Base sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Tbaseline averages per country
# Base data to evaluate functions on (tmx and tmp are Tbaselines)
usdata <- readRDS("tbaseline_countries.rds") %>% filter(country == "USA")
ustmxs = usdata$tmx
ustmps = usdata$tmp
names(ustmxs) <- usdata$Crop
names(ustmps) <- usdata$Crop

# GLOBIOM grid
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

# Functions

lspline_to_betas <- function(fit,cuts) {
  vname = all.vars(formula(fit))[2]
  betas = data.frame(
    sta = -999,
    en = cuts[1],
    intr = fit$coefficients["(Intercept)"],
    slp = fit$coefficients[2]
  )
  # Set a NA slope to zero
  if (is.na(betas$slp[[1]])) {betas$slp[[1]] = 0.0}
  cuts = c(cuts,999) # Add a number to be the last cut
  for (i in 2:(length(cuts))) {
    sta = cuts[i-1]
    en = cuts[i]
    # If the slope is NA, see if the last slope was also NA. 
    # If it was, set to zero, else, set to the last one
    if (is.na(fit$coefficients[i+1])) {
      if (is.na(fit$coefficients[i])) {
        slp = 0.0
      } else {
        # Doing nothing means keep the last value of slp
      }
    } else { # The hopefully most common case, take the number from fit
      slp = fit$coefficients[i+1]
    }
    
    p = data.frame(sta)
    names(p) <- vname
    intr = predict(fit,p) - slp*sta
    # if (is.na(intr) | (slp<0.0)) { # FIXME: We force the last one to be positive
    #   slp = betas$slp[i-1]
    #   intr = betas$intr[i-1]
    # }
    betas = rbind(betas,
                  data.frame(sta=sta,en=en,intr=intr,slp=slp))
    
  }
  # betas = rbind(betas,
  # data.frame(sta=en,en=999,intr=intr,slp=slp))
  rownames(betas) <- NULL
  return(betas)
}
add_lines_betas <- function(betas) {
  for (i in 1:nrow(betas)) {
    pdata = seq(betas$sta[i],betas$en[i],0.1)
    lines(pdata,betas$intr[i] + betas$slp[i]*pdata,col='red',lwd=2.0)
  }
}

# Points shapefile and reference for rasterization
shp = st_read(shpfname)
refrast = raster(reffname)

# Create output folder
dir.create(file.path(outfolder), showWarnings = FALSE, recursive = T)

# Modeled impacts in the US
# moddata = srrefdata
# moddata$dyld <- NA
# moddata$source <- paste0("Av. Tbaseline (",ddversionstring,")")
moddata = data.frame()

# png(filename = paste0(outfolder,"plots.png"), width = 600, height = 300*(length(crops)), unit = "px", pointsize = 16)
par(mfrow = c(length(crops),2))
# crop = crops[1]
for (crop in crops) {
  infname = paste0(infolder,crop,".climdata.rds")
  data = readRDS(infname) %>% filter(year %in% years)
  if (ifiltcnt) {
    data <- filter(data, country %in% filtcnts)
  }
  
  # data["EDD"] <- data[eddvnames[crop]]
  # data["GDD"] <- data[gddvnames[crop]] - data["EDD"]
  if (ieddbound) {
    data["EDD"] <- data[eddvnames[crop]] - data[extvnames[crop]]
  } else {
    data["EDD"] <- data[eddvnames[crop]]
  }
  data["GDD"] <- data[gddvnames[crop]] - data[eddvnames[crop]]
  data["nEDD"] <- data["EDD"]/data["ndays"]
  data["nGDD"] <- data["GDD"]/data["ndays"]
  
  
  # Fit nEDD to Tmax and nGDD to Tavg
  # fitedd = lm(EDD ~ lspline(tmaxmean,cuts), data = data)
  # fitgdd = lm(GDD ~ lspline(tempmean,cuts), data = data)
  fitedd = lm(EDD ~ lspline(tmaxmean,cuts), data = data, weights = areacrop)
  fitgdd = lm(GDD ~ lspline(tempmean,cuts), data = data, weights = areacrop)
  
  # Convert to table of intercepts and slopes using the function lspline_to_betas
  betasedd = lspline_to_betas(fitedd,cuts)
  betasgdd = lspline_to_betas(fitgdd,cuts)
  
  # Add the plots and lines
  with(data,plot(tmaxmean,EDD,main=crop,pch='.', cex = 2.0))
  add_lines_betas(betasedd)
  with(data,plot(tempmean,GDD,main=crop,pch='.', cex = 2.0))
  add_lines_betas(betasgdd)
  
  # Calculate impacts in the US
  dumdata1 =  data.frame(deltat = seq(1,10), tmaxmean = ustmxs[[crop]], tempmean = ustmps[[crop]])
  dumdata2 =  dumdata1
  dumdata2$tmaxmean = dumdata2$tmaxmean + dumdata2$deltat
  dumdata2$tempmean = dumdata2$tempmean + dumdata2$deltat
  
  cropdata = data.frame(deltat = dumdata1$deltat)
  cropdata$dedd = predict(fitedd,newdata = dumdata2) - predict(fitedd,newdata = dumdata1)
  cropdata$dgdd = predict(fitgdd,newdata = dumdata2) - predict(fitgdd,newdata = dumdata1)
  cropdata$dyld = eddsens[[crop]]*cropdata$dedd + gddsens[[crop]]*cropdata$dgdd
  cropdata$dyld = cropdata$dyld*100.0
  cropdata$Crop = crop
  
  moddata = rbind(moddata,cropdata)
  
  # Write each crop and nEDD/nGDD piecewise coefficients and hinges
  # write.csv(betasedd, paste0(outfolder, crop, ".betas.nEDD.csv"), row.names = F)
  # write.csv(betasgdd, paste0(outfolder, crop, ".betas.nGDD.csv"), row.names = F)
}
title(ddversionstring, line = -2, outer = TRUE)
par(mfrow = c(1,1))

# dev.off()

moddata$source <- paste0("Av. Tbaseline (",ddversionstring,")")

compdata = rbind(srrefdata,dplyr::select(moddata,names(srrefdata)))

plotcomp = ggplot(compdata) + 
  geom_line(aes(x=deltat, y = dyld, color = Crop, linetype = source),size = 1.5) +
  ggtitle(ddversionstring)
print(plotcomp)
ggsave(paste0(outfolder,"compare_sr.png"),plot = plotcomp)

# dev.off()



# ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y))
ggplot(data) + geom_point(aes(x = tmaxmean, y = EDD, col = y)) +
  facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90))) +
  scale_color_gradientn(colours = rainbow(5))


ggplot(data) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = trngmean),
             size = 0.5) +
  # facet_wrap(~cut(y,c(-90,-60,-23.15,23.15,30,60,90))) +
  facet_wrap(~cut(y,seq(-90,90,10))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(data) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = y),
             size = 0.5) +
  # facet_wrap(~cut(trngmean,classIntervals(trngmean,n = 4, style="pretty")$brks) + 
               # cut(y,c(-90,-60,-23.15,23.15,30,60,90))) +
  facet_grid(cut(trngmean,c(seq(5,20,5))) ~
               cut(y,seq(-40,50,10))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(filter(data, country %in% c("USA","BRA"))) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = y)) +
  facet_wrap(~country) +
  # facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(filter(data, country %in% c("USA","BRA"))) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = y)) +
  facet_grid(country~cut(y,seq(-60,60,10))) +
  # facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90))) +
  scale_color_gradientn(colours = rainbow(5))




alldata = readRDS(infname)
if (ieddbound) {
  alldata["EDD"] <- alldata[eddvnames[crop]] - alldata[extvnames[crop]]
} else {
  alldata["EDD"] <- alldata[eddvnames[crop]]
}
alldata["GDD"] <- alldata[gddvnames[crop]] - alldata[eddvnames[crop]]
alldata["nEDD"] <- alldata["EDD"]/alldata["ndays"]
alldata["nGDD"] <- alldata["GDD"]/alldata["ndays"]

ggplot(filter(alldata, country %in% c("BRA"))) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = y)) +
  facet_grid(year~cut(y,seq(-60,60,10))) +
  # facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(alldata) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = year),
             size = 0.5) +
  # facet_wrap(~cut(trngmean,classIntervals(trngmean,n = 4, style="pretty")$brks) + 
  # cut(y,c(-90,-60,-23.15,23.15,30,60,90))) +
  facet_grid(cut(trngmean,c(seq(5,20,5))) ~
               cut(y,seq(-40,50,10))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(alldata) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = year),
             size = 0.5) +
  # facet_wrap(~cut(trngmean,classIntervals(trngmean,n = 4, style="pretty")$brks) + 
  # cut(y,c(-90,-60,-23.15,23.15,30,60,90))) +
  facet_grid(zone ~
               cut(y,seq(-40,50,10))) +
  scale_color_gradientn(colours = rainbow(5))

ggplot(alldata) + 
  geom_point(aes(x = tmaxmean, y = EDD, col = year),
             size = 0.5) +
  # facet_wrap(~cut(trngmean,classIntervals(trngmean,n = 4, style="pretty")$brks) + 
  # cut(y,c(-90,-60,-23.15,23.15,30,60,90))) +
  facet_wrap(~zone) +
  scale_color_gradientn(colours = rainbow(5))


# 
# 
# 
# 
# ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y))
# ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y)) +
#   facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90)))
# 
# ggplot(data) + geom_point(aes(x = tmaxmean, y = EDD, col = y))
# plt = ggplot(data) + geom_point(aes(x = tmaxmean, y = EDD, col = y)) +
#   # facet_wrap(~cut(y,seq(-90,90,30)))
#   facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90)))
#   # facet_wrap(~zone)
# plt + scale_color_gradientn(colours = rainbow(5))

#' b0 = "$\\beta_0 + \\frac{1}{2}$"
#' 
#' #'
#' {{b0}}
#' 