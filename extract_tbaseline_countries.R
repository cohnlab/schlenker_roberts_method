require(tidyverse)
require(velox)
require(raster)
require(ncdf4)
require(rasterVis)
require(sf)
require(rgeos)
require(tmap)

outfname = "tbaseline_countries.rds"

#shpfname = 'globiom_grid/worldCR.shp'
# shpfname = 'globiom_grid/COLROW30.shp'
shpfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'

# cropnames = c('Maize','Soybeans','Wheat')
# cropnames = c('Maize','Rice','Soybeans','Wheat','Cotton')

# Mask for CRU (Tbaseline)
crumskfname = "cropmap/cropland2000_grid_cru.nc"
mskval = 0.01

# CRU Tbaselines
crubasefolder  = "../moore_coefficients/CRU_gsavg/"
crutmpfnamesuf = ".gsavg.nc"
crutmxfnamesuf = ".gsavg.nc"

# Open the shapefile in sf
shp = st_read(shpfname)

crop = "Maize"

crutmpfname = paste0(crubasefolder,crop,crutmpfnamesuf)
crutmxfname = paste0(crubasefolder,crop,crutmxfnamesuf)

crumsk = raster(crumskfname, varname = "Band1")
crumsk[crumsk<=mskval] <- NA
crumsk[!is.na(crumsk)] <- 1.0

crutmpnc = raster(crutmpfname, varname = "tmp")
crutmpnc = mask(crutmpnc,crumsk)
crutmxnc = raster(crutmxfname, varname = "tmx")
crutmxnc = mask(crutmxnc,crumsk)

crutmp.vx = velox(crutmpnc)
extcrutmp = crutmp.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
extcrutmp$ID_sp <- shp$ADM0_A3
names(extcrutmp)[names(extcrutmp) == "ID_sp"] <- "country"
names(extcrutmp)[names(extcrutmp) == "out"] <- "tmp"

crutmx.vx = velox(crutmxnc)
extcrutmx = crutmx.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
extcrutmx$ID_sp <- shp$ADM0_A3
names(extcrutmx)[names(extcrutmx) == "ID_sp"] <- "country"
names(extcrutmx)[names(extcrutmx) == "out"] <- "tmx"

outdata = left_join(extcrutmp,extcrutmx,by = "country")

shp %>% left_join(outdata,by = c("ADM0_A3" = "country")) -> shp

subrast = crutmpnc
tit = "tmp"
# breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
breaks = seq(10,40,3)
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plottmp <- tm_shape(subrast) + 
  tm_raster(palette = pal, breaks = breaks, midpoint = NA,
            title = tit) 
plottmp

subrast = crutmxnc
tit = "tmx"
# breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
breaks = seq(10,40,3)
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plottmx <- tm_shape(subrast) + 
  tm_raster(palette = pal, breaks = breaks, midpoint = NA,
            title = tit) 

tm_shape(subrast) + 
  tm_raster(palette = pal, breaks = seq(22,28,1), midpoint = NA,
            title = tit) 

plottmx
tmap_arrange(plottmp,plottmx,nrow=2,asp=2)

#+ 
  # tm_shape(reg) + tm_borders() +
  # tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
  #           legend.outside = T) + 
  # tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
  #           panel.label.size = 1.2,
  #           main.title.position = "center", main.title = tit, main.title.size = 1.0) +
  # tm_facets(nrow = gnrow, ncol = gncol)

saveRDS(outdata, file = outfname)
# write.csv(outdata, file = outfname, row.names = F)

