Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lspline)


calname = "Sacks_ZARC_fill_fill_120d"

infolder  = paste0("sheffield_dfs_vpd/",calname,"/")
outfolder = paste0("vpd_scaling/",calname,"/")

# years = 1991:1993
years = 2001:2008

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

dir.create(file.path(outfolder), showWarnings = FALSE, recursive = TRUE)

# crop = "Soybeans"
for (crop in crops) {
  
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
  
  # Create categories
  cutstemp = seq(15,35,2)
  cutstmax = seq(17,35,2)
  data$tempcat = cut(data$tempmean,cutstemp)
  data$tmaxcat = cut(data$tmaxmean,cutstmax)
  
  
  fit = lm(vpdmean ~ EDD:tempcat,data=data)
  summary(fit)
  data$pred = predict(fit,data)
  
  png(filename = paste0(outfolder,crop,".scatter.png"), width = 1000, height = 800, unit = "px", pointsize = 16)
  
  plt = ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean), size = 0.5) + 
    geom_line(aes(x = EDD, y = pred)) +
    scale_color_gradientn(colours = heat.colors(5, rev = T)) + 
    facet_wrap(~tempcat) +
    ggtitle(crop)
  print(plt)
  dev.off()
  
  png(filename = paste0(outfolder,crop,".barplot.png"), width = 1000, height = 800, unit = "px", pointsize = 16)
  coefs = fit$coefficients[grepl("EDD:",names(fit$coefficients))]
  par(mar=c(5,9,4,2))
  barplot(coefs, horiz = T, las = 2)
  dev.off()
  
  
  outdata = data.frame(cuts = cutstemp, coefs = c(coefs[1],coefs))
  
  outfname = paste0(outfolder, crop, ".vpdscale.csv")
  
  write.csv(outdata, outfname, row.names = F)
  
}

# #TMAX
# fit = lm(vpdmean ~ EDD:tmaxcat,data=data)
# summary(fit)
# data$pred = predict(fit,data)
# ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tmaxmean), size = 0.5) + 
#   geom_line(aes(x = EDD, y = pred)) +
#   scale_color_gradientn(colours = heat.colors(5, rev = T)) + 
#   facet_wrap(~tmaxcat)
# 
# 
# coefs = fit$coefficients[grepl("EDD:",names(fit$coefficients))]
# par(mar=c(5,9,4,2))
# barplot(coefs[-1], horiz = T, las = 2)

# # Tests
# ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean))
# 
# ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) +
#   facet_wrap(~cut(tempmean,seq(15,35,2))) + scale_color_gradientn(colours = rainbow(5))
# 
# ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) +
#   facet_wrap(~zone + cut(tempmean,seq(15,35,5))) + scale_color_gradientn(colours = heat.colors(5, rev = T))
# 
# ggplot(data) + geom_point(aes(x = tempmean, y = vpdmean, col = EDD)) + 
#   scale_color_gradientn(colours = heat.colors(5, rev = T)) +
#   facet_wrap(~zone)
# 
# ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y))
# ggplot(data) + geom_point(aes(x = tempmean, y = GDD, col = y)) +
#   facet_wrap(~cut(y,c(-90,-23.15,23.15,30,90)))
# fit = lm(vpdmean ~ EDD,data=data)
# data$pred = predict(fit,data)
# ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean)) + 
#   geom_line(aes(x = EDD, y = pred)) +
#   scale_color_gradientn(colours = heat.colors(5, rev = T))


fit = lm(vpdmean ~ EDD + EDD:log(tempmean),data=data)
summary(fit)
data$pred = predict(fit,data)

t = seq(15,35)
plot(t, fit$coefficients['EDD'] + fit$coefficients['EDD:log(tempmean)']*log(t))

plt = ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean), size = 0.5) + 
  geom_line(aes(x = EDD, y = pred)) +
  scale_color_gradientn(colours = heat.colors(5, rev = T)) + 
  facet_wrap(~tempcat) +
  ggtitle(crop)


plt = ggplot(data) + geom_point(aes(x = EDD, y = vpdmean, col = tempmean), size = 0.5) + 
  geom_line(aes(x = EDD, y = pred)) +
  scale_color_gradientn(colours = heat.colors(5, rev = T)) + 
  facet_wrap(~tempcat) +
  ggtitle(crop)
print(plt)


