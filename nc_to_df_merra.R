# Just extract the non-NaN points from the computed rasters as a single dataframe per crop

library(tidyverse)
library(raster)
library(sf)

calname = "Sacks_ZARC_fill_fill_120d"
ddversionstring = "merra"

infolder  = paste0("merra_computed//",calname,"/")
outfolder  = paste0("merra_dfs/",calname,"/",ddversionstring,"/")
# infolder  = "xavier_computed/"
# outfolder  = "xavier_dfs/"

# Mask only areas with more agricultural area fraction than mskthresh
mskfname = "cropmap/cropland2000_grid_merra.nc"
mskthresh = 0.01

# Actual area weights per crop
areafpref = "cropmap/crops_merra/"
areafsuf  = ".HarvestedAreaFraction.nc"

zonfname = "GIS/COLROW30_K_G.shp"

dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# crops = c("Soybeans")
# years = 1991:1993

crops = c("Maize","Soybeans","Cotton")
# years = 1991:2008
# years = 2002:2008
years = 2009:2014
# years = 2010:2014

msk = raster(mskfname)
msk[msk<mskthresh] <- NA
msk[msk>mskthresh] <- 1

shpzones = st_read(zonfname)
zonemap = seq(length(levels(shpzones$Class4)))
names(zonemap) <- levels(shpzones$Class4)
shpzones$Class4num = zonemap[shpzones$Class4]
invzonemap = names(zonemap)
names(invzonemap) <- zonemap


for (crop in crops) {
  print(crop)
  allydata = data.frame()
  for (year in years) {
    infname = paste0(infolder,crop,".computed.",year,".nc")
    
    tempmean = raster(infname, varname = "tempmean")
    # Mask just tempmean, na.omit will eliminate other cells in the join phase
    tempmean[is.na(msk)] <- NA
    
    # Rasterize the maps of the zones in the first iteration. Same for 2000 area
    if ((crop == crops[1]) & (year == years[1]) ) {
      rzones = rasterize(shpzones, tempmean, field = "Class4num", fun = "last", background = NA_real_,by = NULL)
      dfzones = as.data.frame(rzones, xy = TRUE)
      dfzones$zone = invzonemap[as.character(dfzones$layer)]
      dfzones$zone[is.na(dfzones$zone)] <- "NODATA"
      dfzones$zone = as.factor(dfzones$zone)
      dfzones$layer <- NULL
      
      areafname = paste0(areafpref,crop,areafsuf)
      areacroprast = raster(areafname)
      areacrop = as.data.frame(areacroprast, xy = TRUE)
      names(areacrop)[3] <- "areacrop"
      areacrop[is.na(areacrop)] <- 0.0
      
    }
    
    tmaxmean = raster(infname, varname = "tmaxmean")
    tminmean = raster(infname, varname = "tminmean")
    # precmean = raster(infname, varname = "precmean")
    # vpdmean = raster(infname, varname = "vpdmean")
    
    trngmean = raster(infname, varname = "trngmean")
    ndays =  raster(infname, varname = "ndays")
    
    tempdist = brick(infname, varname = "tempdist")
    tempgdds = brick(infname, varname = "tempgdds")
    tempdistzvals = getZ(tempdist)
    tempgddszvals = getZ(tempgdds)
    
    tempdist = as.data.frame(tempdist, xy = TRUE)
    colnames(tempdist)[-(1:2)] <- gsub("-","m",paste0("tdi",tempdistzvals))
    
    tempgdds = as.data.frame(tempgdds, xy = TRUE)
    colnames(tempgdds)[-(1:2)] <- gsub("-","m",paste0("edd",tempgddszvals))
    
    tempmean = as.data.frame(tempmean, xy = TRUE)
    tmaxmean = as.data.frame(tmaxmean, xy = TRUE) 
    tminmean = as.data.frame(tminmean, xy = TRUE) 
    # precmean = as.data.frame(precmean, xy = TRUE)
    # vpdmean = as.data.frame(vpdmean, xy = TRUE)
    
    trngmean = as.data.frame(trngmean, xy = TRUE) 
    ndays = as.data.frame(ndays, xy = TRUE) 
    

    
    ydata <- tempmean %>% 
      mutate(pid = group_indices(.,x,y)) %>%
      left_join(tmaxmean, by = c("x","y")) %>% 
      left_join(tminmean, by = c("x","y")) %>% 
      # left_join(precmean, by = c("x","y")) %>%
      # left_join(vpdmean, by = c("x","y")) %>%
      left_join(trngmean, by = c("x","y")) %>% 
      left_join(ndays, by = c("x","y")) %>% 
      left_join(tempdist, by = c("x","y")) %>% 
      left_join(tempgdds, by = c("x","y")) %>%
      left_join(dfzones, by = c("x","y")) %>%
      left_join(areacrop, by = c("x","y")) %>%
      na.omit()
    ydata$year = year
    
    # Actual gdds
    ydata$gdd1030 = ydata$edd10 - ydata$edd30
    
    # Reorder the columns to get the most informative ones first
    firstvars = c("pid","x","y","year")
    varorder = c(firstvars,colnames(ydata)[!(colnames(ydata) %in% firstvars)])
    ydata = ydata[,varorder]
    
    allydata = rbind(allydata,ydata)
  }
  
  outfname = paste0(outfolder,crop,".climdata.rds")
  saveRDS(allydata,file = outfname)
}

