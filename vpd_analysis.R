Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lspline)


calname = "Sacks_ZARC_fill_fill_120d"

infolder  = paste0("sheffield_dfs_vpd/",calname,"/")
outfolder = paste0("gdd_betas/",calname,"/")

# years = 1991:1993
years = 2001:2003

# crops = c("Soybeans")
crops = c("Maize","Soybeans","Cotton")

# Each crop will use different tresholds for GDD and EDD.
# Give lower limit EDD names for GDD, and later GDD can be calculated as EDDl-EDDu
# Values taken from Schlenker and Roberts can be found in digitized/traced.ods
extvnames = c("edd38","edd38","edd38")
eddvnames = c("edd29","edd30","edd32")
gddvnames = c("edd10","edd10","edd15")
names(extvnames) <- crops
names(eddvnames) <- crops
names(gddvnames) <- crops


crop = "Soybeans"

infname = paste0(infolder,crop,".climdata.rds")
data = readRDS(infname) %>% filter(year %in% years)

# data["EDD"] <- data[eddvnames[crop]]
# data["GDD"] <- data[gddvnames[crop]] - data["EDD"]
data["EDD"] <- data[eddvnames[crop]] - data[extvnames[crop]]
data["GDD"] <- data[gddvnames[crop]] - data[eddvnames[crop]]
data["nEDD"] <- data["EDD"]/data["ndays"]
data["nGDD"] <- data["GDD"]/data["ndays"]

# Filter out some zones
data <- data %>% filter(!(zone %in% c("polar","boreal","ocean")))

ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean))

ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) +
  facet_wrap(~cut(tempmean,seq(15,35,2))) + scale_color_gradientn(colours = rainbow(5))

ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) +
  facet_wrap(~zone + cut(tempmean,seq(15,35,5))) + scale_color_gradientn(colours = heat.colors(5, rev = T))

ggplot(data) + geom_point(aes(x = tempmean, y = vpdmean, col = EDD)) + 
  scale_color_gradientn(colours = heat.colors(5, rev = T)) +
  facet_wrap(~zone)

ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y))
ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y)) +
  facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90)))

fit = lm(vpdmean ~ EDD,data=data)
data$pred = predict(fit,data)
ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) + 
  geom_line(aes(x = EDD, y = pred)) +
  scale_color_gradientn(colours = heat.colors(5, rev = T))

fit = lm(vpdmean ~ EDD,data=data)
data$pred = predict(fit,data)
ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) + 
  geom_line(aes(x = EDD, y = pred)) +
  scale_color_gradientn(colours = heat.colors(5, rev = T))




