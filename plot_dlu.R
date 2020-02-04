library(tidyverse)
library(raster)
library(sf)
dlufname = "C:/Users/eriad/Downloads/newAGT1_DeltaLUCs.csv"
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

convdata = read.csv(convfname)
convdata$COLROW30 <- as.character(convdata$final.COLROW30)

dludata = read.csv(dlufname)
dludata <- left_join(dludata,convdata,by = c( "ID" = "final.ID"))

shp = st_read(shpfname)
shp$COLROW30 <- as.character(shp$COLROW30)
ref = raster(reffname)

# y = 2005
# 
# ydata <- dludata %>% filter(year == y) 

shp = left_join(shp,dludata)
# yshp = left_join(shp,ydata, by = c("COLROW30" = "final.COLROW30"))
summary(shp)

# yrast <- rasterize(yshp, ref, field = "dluc", fun = mean, background = NA_real_,
                   # by = NULL)
rast <- rasterize(shp, ref, field = c("X2005","X2050"), fun = mean, background = NA_real_,
                   by = NULL)


plot(shp[c("X2050","geometry")], ylim = c(-33,10), xlim = c(-80,-30))
plot(rast, ylim = c(-33,10), xlim = c(-80,-30))
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
