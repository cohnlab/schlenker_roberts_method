library(tidyverse)
library(raster)
library(sf)
dlufname = "C:/Users/eriad/Downloads/DELTA_LUC_AGT1SSP0.CSV"
shpfname = "GIS/COLROW30.shp"
reffname = "calendars_sacks/Cotton.crop.calendar.fill.nc"
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

convdata = read.csv(convfname)
convdata$COLROW30 <- as.character(convdata$final.COLROW30)

dludata = read.csv(dlufname, header = FALSE)
names(dludata) <- c("country","cell","year","dluc")
dludata <- left_join(dludata,convdata,by = c( "cell" = "final.ID"))

shp = st_read(shpfname)
shp$COLROW30 <- as.character(shp$COLROW30)
ref = raster(reffname)

y = 2005

ydata <- dludata %>% filter(year == y) 

yshp = left_join(shp,ydata)
# yshp = left_join(shp,ydata, by = c("COLROW30" = "final.COLROW30"))
summary(yshp)

yrast <- rasterize(yshp, ref, field = "dluc", fun = mean, background = NA_real_,
                   by = NULL)

plot(yshp[c("dluc","geometry")], ylim = c(-33,10), xlim = c(-80,-30))
plot(yrast, ylim = c(-33,10), xlim = c(-80,-30))
# # data <- data %>% 
# data1 = filter(data,year == 2005)
# data2 = filter(data,year == 2050)
# 
# both = left_join(data1,data2,by = "cell")
# head(both)
# both$diff = both$dlu.y-both$dlu.x
# 
# hist(data2$dlu, probability = T)
# hist(both$diff, freq=F)
# 
# library(ggplot2)
# ggplot(data,aes(x=year, y = dlu, color=cell)) + geom_line()
