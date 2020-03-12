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
  
  # Writing a file for each spec
  outspec1 = data.frame(start = c(-999,cutstemp), end = c(cutstemp,999), coefs = c(coefs[1],coefs,coefs[length(coefs)]))
  write.csv(outspec1, paste0(outfolder, crop, ".vpdscale.spec1.csv"), row.names = F)
  
  # Spec 2
  cutstemp = seq(17,35,2)
  data$tempcat = cut(data$tempmean,cutstemp)
  fit = lm(vpdmean ~ EDD:tempcat,data=data)
  coefs = fit$coefficients[grepl("EDD:",names(fit$coefficients))]
  outspec2 = data.frame(start = c(-999,cutstemp), end = c(cutstemp,999), coefs = c(coefs[1],coefs,coefs[length(coefs)]))
  write.csv(outspec2, paste0(outfolder, crop, ".vpdscale.spec2.csv"), row.names = F)
  
  # Spec 3: Linear interaction
  fit = lm(vpdmean ~ EDD + EDD:tempmean,data=data)
  outspec3 = data.frame(t(fit$coefficients))
  write.csv(outspec3, paste0(outfolder, crop, ".vpdscale.spec3.csv"), row.names = F)
  
  # Spec 3: Linear interaction
  fit = lm(vpdmean ~ EDD + EDD:log(tempmean),data=data)
  outspec4 = data.frame(t(fit$coefficients))
  write.csv(outspec4, paste0(outfolder, crop, ".vpdscale.spec4.csv"), row.names = F)
  
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

# Using Rcpp here to avoid figuring out a faster lookup in R
cppFunction('NumericVector cpp_eval_lut(NumericMatrix x, NumericVector tvec) {
  int nrow = x.nrow(), ncol = x.ncol();
  int vsize = tvec.size();
  int p = 0;
  double t;
  double val;
  NumericVector vout(vsize);
  for (int p = 0; p < vsize; ++p) {
  t = tvec(p);
  for (int i = 0; i < nrow; ++i) {
  if (t >= x(i,0) && t<= x(i,1)) {
    vout(p) = x(i,2);
    break;
  }
  }
  }
  return vout;
}')

# Wrapper for cpp_eval_lut. Forces a single-column data.frame into a vector
eval_lut <- function(betas,tvec) {
  if (is.data.frame(tvec)) {
    tvec = tvec[,]
  }
  return(cpp_eval_lut(as.matrix(betas),as.vector(tvec)))
}


mdata = data.frame(temps = seq(15,35,0.1))

# Spec 1: discrete from 15-35
# Create categories
cutstemp = seq(15,35,2)
data$tempcat = cut(data$tempmean,cutstemp)

fit = lm(vpdmean ~ EDD:tempcat,data=data)
summary(fit)

coefs = fit$coefficients[grepl("EDD:",names(fit$coefficients))]

mtable = data.frame(cutstart = c(-999,cutstemp), cutsend = c(cutstemp,999), coefs = c(coefs[1],coefs,coefs[length(coefs)]))

mdata$spec1 = eval_nxdd(mtable,mdata$temps)

# Spec 2: discrete from 17-35
# Create categories
cutstemp = seq(17,35,2)
data$tempcat = cut(data$tempmean,cutstemp)

fit = lm(vpdmean ~ EDD:tempcat,data=data)
summary(fit)

coefs = fit$coefficients[grepl("EDD:",names(fit$coefficients))]

mtable = data.frame(cutstart = c(-999,cutstemp), cutsend = c(cutstemp,999), coefs = c(coefs[1],coefs,coefs[length(coefs)]))

mdata$spec2 = eval_nxdd(mtable,mdata$temps)


# Spec 3: Linear interaction
fit = lm(vpdmean ~ EDD + EDD:tempmean,data=data)
summary(fit)
mdata$spec3 = fit$coefficients['EDD'] + fit$coefficients['EDD:tempmean']*mdata$temps

# Spec 4: Log interaction
fit = lm(vpdmean ~ EDD + EDD:log(tempmean),data=data)
mdata$spec4 = fit$coefficients['EDD'] + fit$coefficients['EDD:log(tempmean)']*log(mdata$temps)

# Apply scaling
sdata = mdata
refind = which.min(abs(mdata$temps - 22.00))
sdata[,-1] = sweep(mdata[,-1],2,t(mdata[refind,-1]),"/")

# Plot everyone
longmdata <- mdata %>% gather(key = "variable", value = "value", -temps)
longsdata <- sdata %>% gather(key = "variable", value = "value", -temps)

# ggplot(longmdata,aes(x = temps, y = value)) +
  # geom_line(aes(color = variable),size = 1.5) + theme_classic()
ggplot(longsdata,aes(x = temps, y = value)) +
  geom_line(aes(color = variable),size = 1.5) +
  labs(y = "Scaling factor",
       x = "Baseline temperature",
       color = "VPD~EDD model",
       title = crop) +
  scale_color_discrete(labels = c("Binned (15-35)","Binned (17-35)","Linear interaction","Log interaction")) +
  theme_classic(base_size = 14)

