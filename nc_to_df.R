# Just extract the non-NaN points from the computed rasters as a single dataframe per crop

library(tidyverse)
library(raster)

infolder  = "xavier_computed/"
outfolder  = "xavier_dfs/"

crops = c("Soybeans")
years = 2000:2002


for (crop in crops) {
  
  allydata = data.frame()
  for (year in years) {
    infname = paste0(infolder,crop,".computed.",year,".nc")
    
    tempmean = raster(infname, varname = "tempmean")
    tmaxmean = raster(infname, varname = "tmaxmean")
    tminmean = raster(infname, varname = "tminmean")
    precmean = raster(infname, varname = "precmean")
    
    tempdist = brick(infname, varname = "tempdist")
    tempgdds = brick(infname, varname = "tempgdds")
    tempdistzvals = getZ(tempdist)
    tempgddszvals = getZ(tempgdds)
    
    tempdist = as.data.frame(tempdist, xy = TRUE)
    colnames(tempdist)[-(1:2)] <- gsub("-","m",paste0("tdi",tempdistzvals))
    
    tempgdds = as.data.frame(tempgdds, xy = TRUE)
    colnames(tempgdds)[-(1:2)] <- gsub("-","m",paste0("gdd",tempgddszvals))
    
    tempmean = as.data.frame(tempmean, xy = TRUE)
    tmaxmean = as.data.frame(tmaxmean, xy = TRUE) 
    tminmean = as.data.frame(tminmean, xy = TRUE) 
    precmean = as.data.frame(precmean, xy = TRUE)
    
    ydata <- tempmean %>% 
      mutate(pid = group_indices(.,x,y)) %>%
      left_join(tmaxmean, by = c("x","y")) %>% 
      left_join(tminmean, by = c("x","y")) %>% 
      left_join(precmean, by = c("x","y")) %>%
      left_join(tempdist, by = c("x","y")) %>% 
      left_join(tempgdds, by = c("x","y")) %>%
      na.omit()
    ydata$year = year
    
    # Reorder the columns to get the most informative ones first
    firstvars = c("pid","x","y","year")
    varorder = c(firstvars,colnames(ydata)[!(colnames(ydata) %in% firstvars)])
    ydata = ydata[,varorder]
    
    allydata = rbind(allydata,ydata)
  }
  
  outfname = paste0(outfolder,crop,".climdata.rds")
  saveRDS(allydata,file = outfname)
}