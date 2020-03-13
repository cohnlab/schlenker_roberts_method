Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# Function that applies the DD~T curves
source("eval_nxdd.R")

crops = c("Maize","Soybeans","Cotton")

# Base sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

betasfolder = paste0("gdd_betas/Sacks_ZARC_fill_fill_120d/")

caldata = data.frame(tmaxbase = c(27,27,30), deltat = c(1,6,1), deltay = c(-6.4,-55.9,-1))

crop = "Maize"

# Read the crop-specific coefficients for the T-GDD functions
betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))

caldata$deltagdd = eval_nxdd(betasgdd,caldata$tmaxbase + caldata$deltat) - eval_nxdd(betasgdd,caldata$tmaxbase)
caldata$deltaedd = eval_nxdd(betasedd,caldata$tmaxbase + caldata$deltat) - eval_nxdd(betasedd,caldata$tmaxbase)


coef0 = 0.0
coef1 = 1.0


