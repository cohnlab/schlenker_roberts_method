Packages <- c("ggridges","Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

infname = "deploy_tests/results_both_NEW.rds"
outfname = "deploy_tests/Plots_dlu_RCP45_25_50_NEW.pdf"

reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
shpfname = "GIS/COLROW30.shp"

bbox = c(-130,160,-50,65)

crops = c("Soybeans","Maize","Rice","Wheat")

ref = raster(reffname)
shp = st_read(shpfname)
shp = shp %>% dplyr::rename(ID = "COLROW30")

pdf(file = outfname, width = 15, height = 10)

# Loop crops right from the start
# crop = "Soybeans"
for (crop in crops) {
  print(crop)
  alldata = readRDS(infname) %>% filter(Crop == crop) %>%
    # filter(dluNL %in% c(0.25,0.75) & dluL %in% c(0.25,0.75))
    filter(dluNL %in% c(0.25,0.50) & dluL %in% c(0.25,0.50))
  
  
  
  rcpdata = alldata %>% filter(scen == "ssp245")
  
  
  # Make a combination of factors
  
  levels(rcpdata$dluL) <- paste0("L.",levels(rcpdata$dluL))
  levels(rcpdata$dluNL) <- paste0("NL.",levels(rcpdata$dluNL))
  rcpdata$dluscen = interaction(rcpdata$dluL,rcpdata$dluNL)
  rcpdata$dluscen = droplevels(rcpdata$dluscen)
  
  dluscens = levels(rcpdata$dluscen)
  
  
  # Plot for a method
  for (methodstr in c("MR","DD")) {
    print(methodstr)
    
    widercpdata = rcpdata %>% 
      filter(method == methodstr) %>% 
      dplyr::select(-dluL,-dluNL) %>%
      spread(dluscen,dyld)
    
    rcpshp = left_join(shp,widercpdata, by = c("ID"))
    
    rast = rasterize(rcpshp, ref, field = dluscens, fun = mean, background = NA_real_,
                     by = NULL)
    
    breaks = seq(-100,100,20)
    pal = brewer.pal(n = length(breaks), name = "RdBu")
    tit = paste0("RCP4.5 (2050) | ",crop," (",methodstr,")")
    plotobj <- tm_shape(rast, bbox = bbox) + 
      tm_raster(palette = pal, breaks = breaks, midpoint = 0,
                title = expression(paste(Delta,"Yield (%)"))) + 
      tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
                legend.outside = T) + 
      tm_layout(panel.show = T, panel.labels = names(rast),
                panel.label.size = 1.2,
                main.title.position = "center", main.title = tit, main.title.size = 1.0) +
      tm_facets(nrow = 2, ncol = 2)
    print(plotobj)
    
  }
}

dev.off()
# rcpdata %>%
#   object.size %>%
#   print(units = "Mb")
