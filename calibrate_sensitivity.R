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

# caldata = data.frame(tmx = c(27,27,30), X.Temp = c(1,6,1), deltay = c(-6.4,-55.9,-1))
caldata = read.csv("../moore_coefficients/includedstudies_revisions_tbaseline.csv")
names(caldata)[names(caldata) == "Yield.change...."] <- "deltay"

crop = "Maize"

# Read the crop-specific coefficients for the T-GDD functions
betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))

caldata$deltagdd = eval_nxdd(betasgdd,caldata$tmx + caldata$X.Temp) - eval_nxdd(betasgdd,caldata$tmx)
caldata$deltaedd = eval_nxdd(betasedd,caldata$tmx + caldata$X.Temp) - eval_nxdd(betasedd,caldata$tmx)

caldata$gddsens = gddsens[[crop]]

caldata$deltayfrac = 0.01*caldata$deltay

caldata$eddsens = (caldata$deltayfrac - caldata$gddsens*caldata$deltaedd)/caldata$deltaedd

eddstudies = c("Schlenker and Lobell 2010","Schlenker and Roberts 2009")

caldata$modeltype = "Mechanistic"
caldata$modeltype[caldata$statmodel == 1] = "Statistical"
caldata$modeltype[caldata$Reference %in% eddstudies] = "Statistical EDD"

# Filter out positive impacts and specify the crop
caldata %>% 
  filter(Crop1 == crop) %>%
  filter(eddsens <= 0) ->
  data

statdata = filter(data,modeltype != "Mechanistic")

ggplot(filter(data,modeltype != "Mechanistic"),aes(x = tmx, y = eddsens)) +
  geom_point(aes(color = modeltype),size = 3) + theme_classic()

ggplot(statdata,aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4) + theme_classic()


ggplot(filter(statdata, !(Reference %in% c("Moya et al", "Corobov"))),
       aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4) + theme_classic()

hlinevalue = -0.006435774
ggplot(filter(statdata, (Reference %in% c("Schlenker and Roberts 2009"))),
       aes(x = X.Temp, y = eddsens)) +
  geom_point(aes(color = bootstrap.group),size = 4) + theme_classic() +
  xlab("deltaT") + ylim(c(-0.02,-0.006)) + 
  geom_hline(yintercept = hlinevalue) + 
  geom_text(aes(2,hlinevalue,
                label = paste("SR value = ",sprintf("%f", hlinevalue)), vjust = -1))

# View(filter(statdata, (Reference %in% c("Schlenker and Roberts 2009"))))



coef0 = 0.0
coef1 = 1.0


