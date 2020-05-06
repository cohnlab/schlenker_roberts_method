Packages <- c("ggridges","Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

infname = "deploy_tests/results_both_NEW.rds"
outfname = "deploy_tests/Ridges_deltat_NEW.pdf"

reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
shpfname = "GIS/COLROW30.shp"

bbox = c(-130,160,-50,65)

crops = c("Soybeans","Maize","Rice","Wheat")

# Zones
zonesfname = "GIS/COLROW30_K_G.shp"

# File with the area weights for all crops
wgtfname = "GIS/areas.csv"

# File with equivalencies for conversion from the polygon grid to the points grid
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

# Reading zones
zonesshp = st_read(zonesfname)
zonesshp$zonename <- NA
zonesshp$zonecode <- NA
zonesshp$zonename[zonesshp$Class4 == "warm temperate"] <- "TEMP"
zonesshp$zonename[zonesshp$Class4 == "equatorial"] <- "EQUA"
zonesshp$zonename[zonesshp$Class4 == "arid"] <- "ARID"
zonesshp$zonename[is.na(zonesshp$zonename)] <- "OTHR"
zonesshp$zonecode[zonesshp$Class4 == "warm temperate"] <- 1
zonesshp$zonecode[zonesshp$Class4 == "equatorial"] <- 2
zonesshp$zonecode[zonesshp$Class4 == "arid"] <- 3
zonesshp$zonecode[is.na(zonesshp$zonecode)] <- 4
zonecodes = c(1,2,3,4)
names(zonecodes) = c("TEMP","EQUA","ARID","OTHR")

zonesdata = zonesshp[c("COLROW30","zonecode","zonename")] %>% 
  rename(ID = "COLROW30") %>%
  st_set_geometry(NULL)


# Read fractional area weights
wgtdata = read.csv(wgtfname) %>% dplyr::select(-ID) %>% rename(ID = "COLROW30")
wgtdata = gather(wgtdata,Crop,Area,Wheat:Cotton,factor_key = T)


pdf(file = outfname, width = 10, height = 7.5)

# Loop crops right from the start
# crop = "Soybeans"
for (crop in crops) {
  print(crop)
  rcpdata = readRDS(infname) %>% filter(Crop == crop) %>%
    # filter(dluNL %in% c(0.25,0.75) & dluL %in% c(0.25,0.75))
    filter(dluNL %in% c(0.00) & dluL %in% c(0.00)) %>% 
    filter(scen != "ssp245") %>%
    filter(scen %in% c("deltaT1","deltaT3","deltaT5")) %>%
    left_join(zonesdata,by="ID") %>% 
    left_join(filter(wgtdata, Crop == crop), by = c("ID","Crop"))
  
  
  rcpdata$scenmethod = interaction(rcpdata$method,rcpdata$scen,sep = ".")
  rcpdata$scenmethod = droplevels(rcpdata$scenmethod)
  scenmethods = levels(rcpdata$scenmethod)
  
  # Doing the replications
  repdata <- filter(rcpdata,Area >= 0.001)
  repdata$reps = ceiling(repdata$Area/0.002)
  repdata <- repdata[rep(row.names(repdata), repdata$reps), ]
  
  #This is for getting the area fractions
  fracs <- repdata %>% group_by(Crop, zonename) %>% 
    summarise(frac = sum(reps)) #%>% filter(zonename != "OTHR")
  fracs$frac = 100.0*fracs$frac/sum(fracs$frac)
  fracs$zonenamefrac = paste0(fracs$zonename, "(",sprintf("%.1f",fracs$frac),"%)")
  repdata = left_join(fracs,repdata,by = "zonename")
  
  plt <- ggplot(repdata,
                aes(x = dyld, y = method)) + 
    facet_grid(scen ~ zonenamefrac) +
    geom_density_ridges(aes(fill = zonename, alpha = 0.8), 
                        show.legend = F) +
    theme_ridges() + 
    labs(x = "Yield change (%)",
         # y = "Fixed background warming",
         title = paste0(crop," (Yield)"),
         subtitle = "Fixed background warming"
    )
  print(plt)
  
  # plt <- ggplot(repdata,
  #               aes(x = dtemp, y = method)) + 
  #   facet_grid(scen ~ zonenamefrac) +
  #   geom_density_ridges(aes(fill = zonename, alpha = 0.8), 
  #                       show.legend = F) +
  #   theme_ridges() + 
  #   labs(x = "Tmean change (degC)",
  #        y = "Land use change (pp of whole cell)",
  #        title = paste0(crop," (DeltaTmean)"),
  #        subtitle = "RCP4.5"
  #   )
  # print(plt)
  
}

dev.off()


