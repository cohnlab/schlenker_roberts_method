# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

mainfoldermo = "../AgroServYield_coefficients/estimated_global/Feb7/"
mainfoldersr = "rasters_sr_memo_global/Sacks_Zarc_fill_fill_120d/"

yearmo = 2050
yearsr = 2000

vprefmo = "dlu"
vprefsr = "dyppDLU"

dlus = c(0.0,0.25,0.5,1.0)
dluvnamesmo = paste0(vprefmo,sprintf("%.2f",dlus))
dluvnamessr = paste0(vprefsr,sprintf("%.2f",dlus))
names(dluvnamesmo) <- dlus
names(dluvnamessr) <- dlus

# scenmo = "deltaT2degC"
# scensr = "deltaT2"
scenmo = "deltaTssp0"
scensr = "deltaT0"

# Mask, value is set depending on whether its global below
mskfname = "cropmap/cropland2000_grid_cru.nc"
imask = TRUE

iglobal = TRUE
if (iglobal) {
  regfname = "GIS/bbox_world.shp"
  bbox = extent(c(-130,160,-50,65))
  gnrow = length(dlus)
  gncol = 1
  
  msklim = 0.001
} else {
  regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"
  msklim = 0.000
  gnrow = length(dlus)
  gncol = 1
}

brmskfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"

ioverbreaks = TRUE
overbreaksmo = seq(-24,24,3)
overbreaksdif = seq(-50,50,5)


# crop = "Soybeans"
crop = "Maize"

srpath = paste0(mainfoldersr,"/sr_",scensr,"_",yearsr,crop)
srrast = brick(srpath)

if (crop != "Soybeans") {
  mopath = paste0(mainfoldermo,scenmo,"/estimate_range","_",scenmo,"_",yearmo,crop)
  morast = brick(mopath)
} else {
  mopath = paste0(mainfoldermo,scenmo,"/estimate_range","_",scenmo,"_",yearmo,"Rice")
  morastrice = brick(mopath)
  mopath = paste0(mainfoldermo,scenmo,"/estimate_range","_",scenmo,"_",yearmo,"Wheat")
  morastwheat = brick(mopath)
  morast = (morastrice+morastwheat)/2
}


mosub = subset(morast,dluvnamesmo)
srsub = subset(srrast,dluvnamessr)

# Set the names as in Moore
dif = srsub - mosub
names(dif) <- names(mosub)

tit = paste0("SR-Moore\n",crop," | ", scenmo, " | Min:", sprintf("%.2f",min(values(dif), na.rm = T)))
breaks = overbreaksdif
plotdif <- tm_shape(dif, bbox = bbox) + 
  tm_raster(palette = "RdBu", breaks = breaks,
            title = expression(paste(Delta,"Yield (%)"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) + 
  tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
            panel.label.size = 1.2,
            main.title.position = "center", main.title = tit, main.title.size = 1.0) +
  tm_facets(nrow = gnrow, ncol = gncol)

tit = paste0("Moore\n",crop," | ", scenmo, " | Min:", sprintf("%.2f",min(values(mosub), na.rm = T)))
breaks = overbreaksmo
plotmo <- tm_shape(mosub, bbox = bbox) + 
  tm_raster(palette = "RdBu", breaks = breaks,
            title = expression(paste(Delta,"Yield (%)"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) + 
  tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
            panel.label.size = 1.2,
            main.title.position = "center", main.title = tit, main.title.size = 1.0) +
  tm_facets(nrow = gnrow, ncol = gncol)
tmap_arrange(plotmo,plotdif)


