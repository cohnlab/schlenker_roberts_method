# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# File with equivalencies for conversion from the polygon grid to the points grid
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

tempbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmn.csv"
tmaxbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmx.csv"

# Points shapefile, reference raster and region of interest shapefile
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

crops = c("Maize","Soybeans","Cotton")

# AGT file
agtfname = "../AgroServYield_coefficients/Resulting_tables/Feb7/agroservT.csv"

# deltaT scenario
scen = "deltaTssp245"
scenstring = "RCP4.5"
year = 2050
dttfname = paste0("scenario_tables/",scen,".csv")

# Read deltaT scenario just for the year of interest
dttdata <- read.csv(dttfname) %>% filter(ScenYear == year) %>% spread("PARAMETER","Value")

# Read the AGT coefficients
agtdata <- read.csv(agtfname) %>% spread("PARAMETER","Value")

# Read baseline temp data
tempbase <- read.csv(tempbasefname) %>% dplyr::select(-X)
tmaxbase = read.csv(tmaxbasefname) %>% dplyr::select(-X)

# Make it long in respect to crops and join them
tempbase = gather(tempbase,"Crop","tempbase",-ID)
tmaxbase = gather(tmaxbase,"Crop","tmaxbase",-ID)

alldata <- left_join(tempbase,tmaxbase,by=c("ID","Crop")) %>%
  left_join(dttdata, by = c("ID","Crop")) %>%
  left_join(agtdata, by = c("ID")) 
alldata$Crop = as.factor(alldata$Crop)

# Calculate AGT 100% impacts
alldata$dtempfinal = alldata$gccdtemp + alldata$agroservT
alldata$dtmaxfinal = alldata$gccdtmax + alldata$agroservT*2

# Read the shapefile, join it with LUC values and rasterize them
shp = st_read(shpfname)
shp$COLROW30 <- as.character(shp$COLROW30) # Improves compatibility

# Reed the reference raster
ref = raster(reffname)

crop = "Maize"

# Filter baselines for the crop of interest
ydata <-alldata %>% filter(Crop ==  crop)

# Join with the shapefile
yshp = left_join(shp,ydata, by = c("COLROW30" = "ID"))

valid = c("tempbase","tmaxbase","gccdtemp","gccdtmax","dtempfinal","dtmaxfinal") 
yrast <- rasterize(yshp, ref, field = valid, fun = mean, background = NA_real_,
                   by = NULL)



subrast = subset(yrast,"gccdtmax")
tit = paste0(scenstring," ", year, " 0% LU")
#breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
breaks = seq(0,8,0.5)
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plot_gccdtmax <- tm_shape(subrast) +
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Tmax"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) +
  tm_layout(main.title = tit, main.title.size = 1)
plot_gccdtmax

subrast = subset(yrast,"gccdtemp")
tit = paste0(scenstring," ", year, " 0% LU")
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plot_gccdtemp <- tm_shape(subrast) + 
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Tavg"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) +
  tm_layout(main.title = tit, main.title.size = 1)
plot_gccdtemp

subrast = subset(yrast,"dtmaxfinal")
tit = paste0(scenstring," ", year, " 100% LU")
#breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
# breaks = seq(0,8,0.5)
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plot_dtmaxfinal <- tm_shape(subrast) +
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Tmax"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) +
  tm_layout(main.title = tit, main.title.size = 1)
plot_dtmaxfinal

subrast = subset(yrast,"dtempfinal")
tit = paste0(scenstring," ", year, " 100% LU")
pal = rev(brewer.pal(n = length(breaks), name = "Spectral"))
plot_dtempfinal <- tm_shape(subrast) + 
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Tavg"))) + 
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) +
  tm_layout(main.title = tit, main.title.size = 1)
plot_dtempfinal



tmap_arrange(plot_gccdtemp,plot_gccdtmax,plot_dtempfinal,plot_dtmaxfinal, nrow = 2)

