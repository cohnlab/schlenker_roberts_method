# Just extract the non-NaN points from the computed rasters as a single dataframe per crop

library(tidyverse)
library(raster)
library(sf)
library(velox)

calname = "Sacks_ZARC_fill_fill_120d"
ddversionstring = "bcagcfsr"

infolder  = paste0("bcagcfsr_computed//",calname,"/")
outfolder  = paste0("bcagcfsr_dfs/",calname,"/",ddversionstring,"/")
# infolder  = "xavier_computed/"
# outfolder  = "xavier_dfs/"

dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# crops = c("Soybeans")
# years = 2002:2008

crops = c("Maize","Soybeans","Cotton")
# years = 2002:2008
years = 2002:2008
deltats = 0:5 # Must get zero

# Points shapefile
shpfname = 'GIS/COLROW30.shp'

shp = st_read(shpfname)

extract_raster <- function(infname, shp, vname) {
  inrast.vx = velox(raster(infname, varname = vname))
  # ext = inrast.vx$extract_points(sp = shp, df = TRUE, small = TRUE, fun = function(x) mean(x, na.rm = TRUE))
  ext = inrast.vx$extract_points(sp = shp)
  return(ext[,1])
}

for (crop in crops) {
  print(crop)
  allydata = data.frame()
  for (year in years) {
    for (deltat in deltats) {
      infname = paste0(infolder,crop,".computed.deltat.",deltat,".",year,".nc")
      
      dtdata = st_set_geometry(shp,NULL)
      dtdata$year = year
      dtdata$deltat = deltat
      
      dtdata$tempmean = extract_raster(infname,shp,"tempmean")
      dtdata$tmaxmean = extract_raster(infname,shp,"tmaxmean")
      dtdata$tminmean = extract_raster(infname,shp,"tminmean")
      dtdata$trngmean = extract_raster(infname,shp,"trngmean")
      dtdata$ndays = extract_raster(infname,shp,"ndays")
      
      inbrick = brick(infname, varname = "tempgdds")
      inbrick.vx = velox(inbrick)
      
      ext = as.data.frame(inbrick.vx$extract_points(sp = shp))
      names(ext) <- paste0("edd",getZ(inbrick))
      dtdata = cbind(dtdata,ext)
      
      
      allydata = rbind(allydata,dtdata)
    }
  }
  
  outfname = paste0(outfolder,crop,".climdata.deltat.rds")
  saveRDS(allydata,file = outfname)
}

