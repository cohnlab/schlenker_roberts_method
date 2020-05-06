Packages <- c("ggridges","Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

infname = "deploy_tests/results_both_NEW.rds"
outfname = "deploy_tests/Plots_dlus_NEW.pdf"

reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
shpfname = "GIS/COLROW30.shp"

bbox = c(-130,160,-50,65)

crops = c("Soybeans","Maize","Rice","Wheat")

ref = raster(reffname)
shp = st_read(shpfname)
shp = shp %>% dplyr::rename(ID = "COLROW30")

pdf(file = outfname, width = 10, height = 7.5)

# Loop crops right from the start
# crop = "Soybeans"
for (crop in crops) {
  print(crop)
  rcpdata = readRDS(infname) %>% filter(Crop == crop) %>%
    filter(dluNL %in% c(0.00,0.50) & dluL %in% c(0.00,0.50)) %>%
    # filter(dluNL %in% c(0.00) & dluL %in% c(0.00)) %>% 
    filter(scen == "ssp245") #%>%
    # filter(scen %in% c("deltaT1","deltaT3","deltaT5"))
  
  levels(rcpdata$dluL) <- paste0("L.",levels(rcpdata$dluL))
  levels(rcpdata$dluNL) <- paste0("NL.",levels(rcpdata$dluNL))
  
  rcpdata$scenmethod = interaction(rcpdata$method,rcpdata$dluL,rcpdata$dluNL,sep = ".")
  rcpdata$scenmethod = droplevels(rcpdata$scenmethod)
  scenmethods = levels(rcpdata$scenmethod)
  
  # rcpdata %>%
    # dplyr::select(-method,-dtemp,-dluL,-dluNL) %>% summary
  
  widercpdata = rcpdata %>%
    dplyr::select(-method,-scen,-dtemp,-dluL,-dluNL) %>%
    spread(scenmethod,dyld) 
  
  rcpshp = left_join(shp,widercpdata, by = c("ID"))
  
  rast = rasterize(rcpshp, ref, field = scenmethods, fun = mean, background = NA_real_,
                   by = NULL)
  
  breaks = seq(-100,100,20)
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  tit = paste0("RCP4.5 (2050) | ",crop)
  # tit = crop
  plotobj <- tm_shape(rast, bbox = bbox) + 
    tm_raster(palette = pal, breaks = breaks, midpoint = 0,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = names(rast),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = 4, ncol = 2)
  print(plotobj)
  
  for (methodstr in c("DD","MR")) {
    subrast = raster::subset(rast, str_subset(names(rast),methodstr))
    tit = paste0(methodstr," | RCP4.5 (2050) | ",crop)
    # tit = crop
    plotobj <- tm_shape(subrast, bbox = bbox) + 
      tm_raster(palette = pal, breaks = breaks, midpoint = 0,
                title = expression(paste(Delta,"Yield (%)"))) + 
      tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
                legend.outside = T) + 
      tm_layout(panel.show = T, panel.labels = names(subrast),
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
