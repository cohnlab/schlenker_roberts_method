Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lspline)


calname = "Sacks_ZARC_fill_fill_120d"
# ddversionstring = "merra"
# ddversionstring = "wgt_crop_unbound"
# ddversionstring = "dyncuts_agcfsr_unbound_usa"
ddversionstring = "agcfsr_unbound_usa"

infolder  = paste0("agcfsr_dfs/",calname,"/","agcfsr","/")
# infolder  = paste0("sheffield_dfs_vpd/",calname,"/",ddversionstring,"/")
# infolder  = paste0("merra_dfs/",calname,"/",ddversionstring,"/")
outfolder = paste0("gdd_betas/",calname,"/",ddversionstring,"/")

# Filter specific countries?
ifiltcnt = TRUE
filtcnts = c("USA")

# years = 2009:2014
years = 2002:2008
# years = 2002:2003
# years = 1991:2008

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
cutstemp = seq(9,37,2)
cutstmax = seq(9,37,2)



# Countries shapefile for filtering
shpfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'

shp = st_read(shpfname)[c("ADM0_A3","geometry")]


# Functions

# FIXME: This function DOES allow for a decreasing last piece,
# problems with rank-defficiency might arise
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
filter_cuts <- function(data,incuts,vname) {  
  nmin = 100 # Minimum points to estimate the first and last pieces
  incuts = cutstemp
  outcuts = incuts
  vname = "tempmean"
  # Forward pass
  n = 0
  c = 0
  while(n <= nmin) {
    c = c+1
    cdata = filter(data,data[[vname]] < incuts[c])
    n = nrow(cdata)
  }
  outcuts = incuts[c:length(incuts)]
  # Backward pass
  n = 0
  c = length(outcuts)+1
  while(n <= nmin) {
    c = c-1
    cdata = filter(data,data[[vname]] > incuts[c])
    n = nrow(cdata)
  }
  outcuts = outcuts[1:c]
  return(outcuts)
}

add_lines_betas <- function(betas) {
  for (i in 1:nrow(betas)) {
    pdata = seq(betas$sta[i],betas$en[i],0.1)
    lines(pdata,betas$intr[i] + betas$slp[i]*pdata,col='red',lwd=2.0)
  }
}


# Create output folder
dir.create(file.path(outfolder), showWarnings = FALSE, recursive = T)

png(filename = paste0(outfolder,"plots.png"), width = 600, height = 300*(length(crops)), unit = "px", pointsize = 16)
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
  fitedd = lm(EDD ~ lspline(tmaxmean,cutstmax), data = data, weights = areacrop)
  fitgdd = lm(GDD ~ lspline(tempmean,cutstemp), data = data, weights = areacrop)
  
  # Convert to table of intercepts and slopes using the function lspline_to_betas
  betasedd = lspline_to_betas(fitedd,cutstmax)
  betasgdd = lspline_to_betas(fitgdd,cutstemp)
  
  # Add the plots and lines
  with(data,plot(tmaxmean,EDD,main=crop,pch='.', cex = 2.0))
  points(data$tmaxmean,predict(fitedd,data), col = "blue", pch = '.',cex = 4.0)
  add_lines_betas(betasedd)
  with(data,plot(tempmean,GDD,main=crop,pch='.', cex = 2.0))
  points(data$tempmean,predict(fitgdd,data), col = "blue", pch = '.',cex = 4.0)
  add_lines_betas(betasgdd)
  

  
  # Write each crop and nEDD/nGDD piecewise coefficients and hinges
  write.csv(betasedd, paste0(outfolder, crop, ".betas.nEDD.csv"), row.names = F)
  write.csv(betasgdd, paste0(outfolder, crop, ".betas.nGDD.csv"), row.names = F)
}
dev.off()

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