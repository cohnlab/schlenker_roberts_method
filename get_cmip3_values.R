require(tidyverse)
require(velox)
require(raster)
require(ncdf4)
require(rasterVis)
require(sf)
require(rgeos)

# infolder = "input_cmip_deltat_ens/"
outfolder = "test_output_tables/"

#shpfname = 'globiom_grid/worldCR.shp'
# shpfname = 'globiom_grid/COLROW30.shp'
# shpfname = 'GIS/Countries_WGS84.shp'
shpfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'

# cropnames = c('Maize','Soybeans','Wheat')
# cropnames = c('Maize','Rice','Soybeans','Wheat','Cotton')

# Mask for cmip3
mskfname = "cropmap/cropland2000_grid_cmip3_flip.nc"
mskval = 0.01

# Mask for CRU (Tbaseline)
crumskfname = "cropmap/cropland2000_grid_cru.nc"

# CRU Tbaselines
crutmpfname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
crutmxfname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmx.nc"


# NOT IMPLEMENTED!
varnames = c("tas")
outvarnames = c("gccdtemp")
# varnames = c("tas","tasmax","tasmin")
# outvarnames = c("gccdtemp","gccdtmax","gccdtmin")
names(outvarnames) <- varnames



# rcps = c('ssp245','ssp585')

# Create output directory
dir.create(outfolder, showWarnings = FALSE)

# Open the shapefile in sf
shp = st_read(shpfname)

# crop = 'Maize'
# rcp = 'ssp245'


fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
} 

crop = "Soybeans"

# infname = paste0(crop,'.',rcp,'.delta.ensmean.nc')
# scen = 'shortA2'
scen = 'sresa2'
infname = paste0("data/test/",scen,"_hist_flip.nc")
print(infname)

msk = raster(mskfname, varname = "Band1")
msk[msk<=mskval] <- NA
msk[!is.na(msk)] <- 1.0

# Open the NetCDF file (year,lat,lon) as a Rasterbrick
nc = brick(infname)
nc = mask(nc,msk)

# Get a VeloxRaster version of it
nc.vx = velox(nc)

# Extract all years at once and get the dates as variable names from the original RasterBrick
#ext = nc.vx$extract(sp = shp, fun = function(x) mean(x, na.rm = TRUE), small = TRUE)
# ext = nc.vx$extract_points(sp = shp)
ext = nc.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
# ext$ID_sp <- shp$CNTRY_NAME
ext$ID_sp <- shp$ADM0_A3
names(ext)[names(ext) == "ID_sp"] <- "country"
colnames(ext)[-1] <- names(nc)

# Convert to long format and extract month and year
data = gather(ext,key = "date", value = "tas", -country)
data$date = as.Date(data$date, format = "X%Y.%m.%d")
data$month = as.numeric(format(data$date,"%m"))
data$year = as.numeric(format(data$date,"%Y"))


# Build the output file name

outfname = paste0(outfolder,'/',scen,'.all.rds')

saveRDS(data,outfname)

# Do the same for Tbaseline
crumsk = raster(crumskfname, varname = "Band1")
crumsk[crumsk<=mskval] <- NA
crumsk[!is.na(crumsk)] <- 1.0

crutmpnc = brick(crutmpfname, varname = "tmp")
crutmpnc = mask(crutmpnc,crumsk)
crutmxnc = brick(crutmxfname, varname = "tmx")
crutmxnc = mask(crutmxnc,crumsk)

crutmp.vx = velox(crutmpnc)
extcrutmp = crutmp.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
extcrutmp$ID_sp <- shp$CNTRY_NAME
names(extcrutmp)[names(extcrutmp) == "ID_sp"] <- "country"
colnames(extcrutmp)[-1] <- names(crutmpnc)
datacrutmp = gather(extcrutmp,key = "date", value = "tas", -country)
datacrutmp$date = as.Date(datacrutmp$date, format = "X%Y.%m.%d")
datacrutmp$month = as.numeric(format(datacrutmp$date,"%m"))

crutmx.vx = velox(crutmxnc)
extcrutmx = crutmx.vx$extract(sp = shp, fun = function(x) mean(x,na.rm=T), df = TRUE,small = TRUE)
extcrutmx$ID_sp <- shp$CNTRY_NAME
names(extcrutmx)[names(extcrutmx) == "ID_sp"] <- "country"
colnames(extcrutmx)[-1] <- names(crutmxnc)
datacrutmx = gather(extcrutmx,key = "date", value = "tas", -country)
datacrutmx$date = as.Date(datacrutmx$date, format = "X%Y.%m.%d")
datacrutmx$month = as.numeric(format(datacrutmx$date,"%m"))

usemon = c(5,6,7,8)

datacrutmp %>% 
  filter(country == "United States") %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) %>% 
  dplyr::select(value) %>% as.numeric

datacrutmx %>% 
  filter(country == "United States") %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) %>% 
  dplyr::select(value) %>% as.numeric


usemon = c(11,12,1,2,3,4,5)
safcountries = c("South Africa","Mozambique","Tanzania, United Republic of","Zimbabwe")

datacrutmp %>% 
  filter(country %in% safcountries) %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) 

datacrutmx %>% 
  filter(country %in% safcountries) %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) 

usemon = c(10,11,12,1,2,3,4,5)
datacrutmp %>% 
  filter(country == "Brazil") %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) 

datacrutmx %>% 
  filter(country  == "Brazil") %>%
  filter(month %in% usemon) %>%
  group_by(country) %>% summarize(value = mean(tas,na.rm = T)) 



# datacrutmx %>% 
#   filter(country == "United States") %>%
#   filter(month %in% usemon) %>%
#   group_by(country) %>% summarize(value = mean(tas,na.rm = T)) %>% 
#   select(value) %>% as.numeric



# Get Tbaseline



# bbox = extent(c(230,300,25,50))
# us = crop(nc,bbox)
# plot(us, main = paste0(scen," deltaT"))
# plot(nc, main = "poi")

# df %>% filter(is.na(X2015))
# 
# for (n in colnames(df)) {
#   shp[n] = df[n]
# }
# # file_name <- system.file("C:/Users/eriad/Downloads/mesoregiao/mesoregiao.shp", package="sf")
# ext = ano$extract(sp = shp, fun = function(x) mean(x, na.rm = TRUE), small = TRUE)
# # ext = nc.vx$extract(sp = shp, fun = function(x) mean(x, na.rm = TRUE))
# head(ext)
# shp['lala'] = ext
# 
# # lala = select(shp[drop=TRUE],-c(geometry))
# 
# st_write(shp,'globiom_grid/test.shp')


data %>% filter(country == "USA", month %in% c(5,6,7,8)) -> usdata
