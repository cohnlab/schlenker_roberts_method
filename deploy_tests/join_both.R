Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# mrfname = "deploy_tests/MRdluresults.rds"\
mrfname = "deploy_tests/results_table_MR.rds"
ddfname = "deploy_tests/results_table_DD_NEW.rds"

outfname = "deploy_tests/results_both_NEW.rds"

# Drop dtmax if it's still there
mrdata = readRDS(mrfname) %>% dplyr::select(-dtmax)
dddata = readRDS(ddfname) %>% dplyr::select(-dtmax)
# mrdata = readRDS(mrfname)
# dddata = readRDS(ddfname)

mrdata$method = "MR"
dddata$method = "DD"

# Build Soybeans with Moore estimates
tempric = filter(mrdata,Crop == "Rice")
names(tempric)[names(tempric) == "dyld"] <- "V1"
names(tempric)[names(tempric) == "dtemp"] <- "dtemp1"
# names(tempric)[names(tempric) == "dtmax"] <- "dtmax1"
tempric$Crop <- NULL
tempwht = filter(mrdata,Crop == "Wheat")
names(tempwht)[names(tempwht) == "dyld"] <- "V2"
names(tempwht)[names(tempwht) == "dtemp"] <- "dtemp2"
# names(tempwht)[names(tempwht) == "dtmax"] <- "dtmax2"
tempwht$Crop <- NULL
temp = left_join(tempric,tempwht, by = c("ID", "scen", "dluL", "dluNL","method"))
rm(tempric,tempwht)
temp$dyld = (temp$V1 + temp$V2)/2.0
temp$dtemp = (temp$dtemp1 + temp$dtemp2)/2.0
# temp$dtmax = (temp$dtmax1 + temp$dtmax2)/2.0
# temp = temp %>% select(-V1,-V2,-dtemp1,-dtemp2,-dtmax1,-dtmax2) 
temp = temp %>% select(-V1,-V2,-dtemp1,-dtemp2) 
  temp$Crop = "Soybeans"
mrdata = rbind(mrdata,temp)

# Build Rice and Wheat with DD estimates
temp = filter(dddata,Crop == "Soybeans")
temp$Crop = "Rice"
# dddata = rbind(dddata,temp)
dddata = bind_rows(dddata,temp)
temp$Crop = "Wheat"
# dddata = rbind(dddata,temp)
dddata = bind_rows(dddata,temp)

# Make sure everyone that should be a factor is a factor
# Do this for both before joining since it can speed up the join
mrdata$scen <- as.factor(mrdata$scen)
mrdata$Crop <- as.factor(mrdata$Crop)
mrdata$dluL <- as.factor(mrdata$dluL)
mrdata$dluNL <- as.factor(mrdata$dluNL)
dddata$scen <- as.factor(dddata$scen)
dddata$Crop <- as.factor(dddata$Crop)
dddata$dluL <- as.factor(dddata$dluL)
dddata$dluNL <- as.factor(dddata$dluNL)

# Join both datasets and remove junk
rm(temp)
bothdata = bind_rows(mrdata,dddata)
rm(mrdata,dddata)


bothdata$method = as.factor(bothdata$method)

# Save 
saveRDS(bothdata,outfname)

