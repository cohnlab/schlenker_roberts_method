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

# crops = c('Maize','Soybeans','Wheat')
crops = c('Maize','Rice','Soybeans','Wheat','Cotton')

# Mask for CRU (Tbaseline)
# crumskfname = "cropmap/cropland2000_grid_cru.nc"
# mskval = 0.01

# Masks for individual crops
cropmskfpref = "cropmap/crops_cru/"
cropmskfsuf = ".HarvestedAreaFraction.nc"
mskval = 0.01 # Not used anymore, weighting now

# CRU Tbaselines
crubasefolder  = "../moore_coefficients/CRU_gsavg/"
crutmpfnamesuf = ".gsavg.nc"
crutmxfnamesuf = ".gsavg.nc"

# Open the shapefile in sf
shp = st_read(shpfname)


outdata = data.frame()
# crop = "Maize"
for (crop in crops) {
  crutmpfname = paste0(crubasefolder,crop,crutmpfnamesuf)
  crutmxfname = paste0(crubasefolder,crop,crutmxfnamesuf)
  
  crutmpnc = raster(crutmpfname, varname = "tmp")
  crutmxnc = raster(crutmxfname, varname = "tmx")
  
  crumskfname = paste0(cropmskfpref,crop,cropmskfsuf)
  
  crumsk = raster(crumskfname, varname = "Band1")
  # crumsk[crumsk<=mskval] <- NA
  # crumsk[!is.na(crumsk)] <- 1.0
  
  crutmpnc = mask(crutmpnc,crumsk)
  crutmxnc = mask(crutmxnc,crumsk)
  
  crumsk.vx = velox(crumsk)
  crumsk.ext = crumsk.vx$extract(sp = shp, df = TRUE,small = TRUE)
  
  crutmp.vx = velox(crutmpnc)
  crutmp.ext = crutmp.vx$extract(sp = shp, df = TRUE,small = TRUE)
  crutmp.ext$wgt = crumsk.ext$do.call..rbind...out.
  crutmp.ext %>% 
    group_by(ID_sp) %>% 
    summarise(out = stats::weighted.mean(do.call..rbind...out.,wgt,na.rm=T)) ->
    extcrutmp
  # extcrutmp = crutmp.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
  extcrutmp$ID_sp <- shp$ADM0_A3
  names(extcrutmp)[names(extcrutmp) == "ID_sp"] <- "country"
  names(extcrutmp)[names(extcrutmp) == "out"] <- "tmp"
  
  crutmx.vx = velox(crutmxnc)
  crutmx.ext = crutmx.vx$extract(sp = shp, df = TRUE,small = TRUE)
  crutmx.ext$wgt = crumsk.ext$do.call..rbind...out.
  crutmx.ext %>% 
    group_by(ID_sp) %>% 
    summarise(out = stats::weighted.mean(do.call..rbind...out.,wgt,na.rm=T)) ->
    extcrutmx
  # extcrutmx = crutmx.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
  extcrutmx$ID_sp <- shp$ADM0_A3
  names(extcrutmx)[names(extcrutmx) == "ID_sp"] <- "country"
  names(extcrutmx)[names(extcrutmx) == "out"] <- "tmx"
  
  cropoutdata = left_join(extcrutmp,extcrutmx,by = "country")
  cropoutdata$Crop = crop
  
  outdata = rbind(outdata,cropoutdata)
}

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

