# Just extract the non-NaN points from the computed rasters as a single dataframe per crop

library(tidyverse)
library(raster)

infolder  = "sheffield_computed/"
outfolder  = "sheffield_dfs/"
# infolder  = "xavier_computed/"
# outfolder  = "xavier_dfs/"

# Mask only areas with more agricultural area fraction than mskthresh
mskfname = "cropmap/cropland2000_grid_sheffield.nc"
mskthresh = 0.01

dir.create(outfolder, showWarnings = FALSE)

# crops = c("Soybeans")
# years = 2000:2002

crops = c("Maize","Soybeans","Rice","Wheat")
years = 1991:2008

msk = raster(mskfname)
msk[msk<mskthresh] <- NA
msk[msk>mskthresh] <- 1

for (crop in crops) {
  print(crop)
  allydata = data.frame()
  for (year in years) {
    infname = paste0(infolder,crop,".computed.",year,".nc")
    
    tempmean = raster(infname, varname = "tempmean")
    # Mask just tempmean, na.omit will eliminate other cells in the join phase
    tempmean[is.na(msk)] <- NA
    
    tmaxmean = raster(infname, varname = "tmaxmean")
    tminmean = raster(infname, varname = "tminmean")
    precmean = raster(infname, varname = "precmean")
    
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
    precmean = as.data.frame(precmean, xy = TRUE)
    
    trngmean = as.data.frame(trngmean, xy = TRUE) 
    ndays = as.data.frame(ndays, xy = TRUE) 
    
    ydata <- tempmean %>% 
      mutate(pid = group_indices(.,x,y)) %>%
      left_join(tmaxmean, by = c("x","y")) %>% 
      left_join(tminmean, by = c("x","y")) %>% 
      left_join(precmean, by = c("x","y")) %>%
      left_join(trngmean, by = c("x","y")) %>% 
      left_join(ndays, by = c("x","y")) %>% 
      left_join(tempdist, by = c("x","y")) %>% 
      left_join(tempgdds, by = c("x","y")) %>%
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