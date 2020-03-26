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

# Base data to evaluate functions on (tmx and tmp are Tbaselines)
basedata = readRDS("tbaseline_countries.rds")

# caldata = data.frame(tmx = c(27,27,30), X.Temp = c(1,6,1), deltay = c(-6.4,-55.9,-1))
caldata = read.csv("../moore_coefficients/includedstudies_revisions_tbaseline.csv")
names(caldata)[names(caldata) == "Yield.change...."] <- "deltay"

crop = "Maize"

# Read the crop-specific coefficients for the T-GDD functions
betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))



caldata$deltagdd = eval_nxdd(betasgdd,caldata$tmp + caldata$X.Temp) - eval_nxdd(betasgdd,caldata$tmp)
caldata$deltaedd = eval_nxdd(betasedd,caldata$tmx + caldata$X.Temp) - eval_nxdd(betasedd,caldata$tmx)

caldata$gddsens = gddsens[[crop]]

caldata$deltayfrac = 0.01*caldata$deltay

caldata$eddsens = (caldata$deltayfrac - caldata$gddsens*caldata$deltagdd)/caldata$deltaedd

eddstudies = c("Schlenker and Lobell 2010","Schlenker and Roberts 2009")

caldata$modeltype = "Mechanistic"
caldata$modeltype[caldata$statmodel == 1] = "Statistical"
caldata$modeltype[caldata$Reference %in% eddstudies] = "Statistical EDD"

ggplot(filter(caldata,X.Temp >= 5),aes(x = tmp, y = deltay)) +
  geom_point(aes(color = modeltype),size = 2) + theme_classic()

ggplot(filter(caldata,X.Temp >= 4),aes(x = tmp, y = deltay)) +
  geom_point(aes(color = X.Temp, shape = modeltype),size = 2) + theme_classic()


ggplot(filter(caldata,statmodel == 1 & X.Temp >= 3),aes(x = tmp, y = deltay)) +
  geom_point(aes(shape = modeltype, color = Reference),size = 2) + theme_classic()

# Filter out positive impacts and specify the crop
caldata %>% 
  filter(Crop1 == crop) %>%
  filter(eddsens <= 0) ->
  data

statdata = filter(data,modeltype != "Mechanistic")

ggplot(filter(data,eddsens >= -0.01),aes(x = tmx, y = eddsens)) +
  geom_point(aes(color = modeltype),size = 4) + theme_classic()


ggplot(filter(data,modeltype == "Mechanistic"),aes(x = tmx, y = eddsens)) +
  geom_point(aes(color = modeltype),size = 4) + theme_classic()

ggplot(filter(data,modeltype != "Mechanistic"),aes(x = tmx, y = eddsens)) +
  geom_point(aes(color = modeltype),size = 3) + theme_classic()

ggplot(statdata,aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4) + theme_classic()


ggplot(filter(statdata, !(Reference %in% c("Moya et al", "Corobov"))),
       aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4) + theme_classic()

ggplot(filter(statdata, tmx >= 30 & (Reference %in% c("Lobell et al. 2008"))),
       aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4)

hlinevalue = -0.006435774
ggplot(filter(caldata, (Reference %in% c("Schlenker and Roberts 2009"))),
       aes(x = X.Temp, y = eddsens)) +
  geom_point(aes(color = bootstrap.group),size = 4) + theme_classic() +
  xlab("deltaT") + #ylim(c(-0.02,-0.006)) + 
  geom_hline(yintercept = hlinevalue) + 
  geom_text(aes(2,hlinevalue,
                label = paste("SR value = ",sprintf("%f", hlinevalue)), vjust = -1))

ggplot(statdata,
       aes(x = tmp, y = X.Temp)) + xlab("Tbaseline") + ylab("deltaT") +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4)


ggplot(filter(statdata, tmx >= 30 & (Reference %in% c("Lobell et al. 2008"))),
       aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4)


ggplot(filter(statdata, X.Temp >= 1 & (Reference %in% c("Lobell et al. 2008"))),
       aes(x = tmx, y = deltay)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4)


# View(filter(statdata, (Reference %in% c("Schlenker and Roberts 2009"))))

# View(filter(statdata, (Reference %in% c("Lobell et al. 2008"))))


# Testing Tmean Tmax relationship
ggplot(basedata) + geom_point(aes(x = tmp, y = tmx))
summary(lm(tmx ~ tmp, basedata))
summary(lm(tmp ~ tmx, basedata))
fittmx = lm(tmp ~ tmx, basedata)

# Synthetic data for plotting
pdata = data.frame(tmx = seq(15,36,3))
pdata$tmp = fittmx$coefficients[[1]] + pdata$tmx*fittmx$coefficients[[2]]
dumdata = pdata
pdata = data.frame()
for (dt in seq(0,6,0.5)) {
  dumdata$deltat = dt
  pdata = rbind(pdata,dumdata)
}
pdata$dgdd = eval_nxdd(betasgdd,pdata$tmp+pdata$deltat) - eval_nxdd(betasgdd,pdata$tmp)
pdata$dedd = eval_nxdd(betasedd,pdata$tmx+pdata$deltat) - eval_nxdd(betasedd,pdata$tmx)

pdata$dyld = pdata$dgdd*gddsens[[crop]] + pdata$dedd*eddsens[[crop]]

ggplot(pdata) + geom_line(aes(x = deltat, y = dgdd, color = as.factor(tmx)))
ggplot(pdata) + geom_line(aes(x = deltat, y = dedd, color = as.factor(tmx)))

ggplot(pdata) + 
  geom_line(aes(x = deltat, y = dyld, color = as.factor(tmx)), size = 1) + ylim(-1,0.1)
  


# coef0 = 0.0
# coef1 = -1.0
# indata = pdata
eval_scaled <- function(indata,coef0,coef1) {
  indata$dgdd = eval_nxdd(betasgdd,indata$tmp+indata$deltat) - eval_nxdd(betasgdd,indata$tmp)
  indata$dedd = eval_nxdd(betasedd,indata$tmx+indata$deltat) - eval_nxdd(betasedd,indata$tmx)
  
  indata$scalesens = coef0 + coef1*log(indata$tmx)
  indata$noscaeldyld = indata$dgdd*gddsens[[crop]] + indata$dedd*eddsens[[crop]]
  indata$dyld = indata$dgdd*gddsens[[crop]] + indata$dedd*eddsens[[crop]]*indata$scalesens
  return(indata$dyld)
}

eval_fixscale_point <- function(tmx,deltat,scalesens) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}

eval_noscale_point <- function(tmx,deltat) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]
  return(dyld)
}

eval_scaled_point <- function(tmx,deltat,coef0,coef1) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  k = 1
  a = 0.39
  m = coef0
  b = coef1
  scalesens = k - (k-a)/(1 + exp(-b*(tmx-m)))
  # scalesens = coef0 + coef1*log(tmx)
  # scalesens = coef0 + coef1*(tmx)
  # scalesens = 1
  noscaeldyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}

eval_scaled_step_point <- function(tmx,deltat) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  scalesens = 4.0832 - 0.1159*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(0.42,scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}

coef0 = 0.0
coef1 = -1.0
pdata$dyld = eval_scaled(pdata,0,1)

eval_scaled_point(26.6,seq(1,10),0,0)

plot(seq(1,10),eval_scaled_point(26.6,seq(1,10),0,0))

plot(seq(1,10),eval_noscale_point(26.6,seq(1,10)))

eval_noscale_point(26.6,6)*1.25
eval_noscale_point(31.0,6)

eval_fixscale_point(31.0,6,0.42)

llim = 0.42
fitlinear = lm(sens~t,data = data.frame(sens = c(1,0.49),t = c(26.6,31.0)))
x = seq(15,40,0.5)
y = fitlinear$coefficients[[1]] + fitlinear$coefficients[[2]]*x
y = pmin(1,y)
y = pmax(llim,y)
plot(x,y)

eval_scaled_point(26.6,6,0,0)
eval_scaled_point(31,6,0,0)

fdata = data.frame(tmx = c(26.6,31.0,32.0), 
                   deltat = c(6,6,6), 
                   dyld = c(-0.54,-0.54*1.25,-0.54*1.25))

scoef0 = 0.0
scoef1 = -1.0
plot(seq(15,40),scoef0 + scoef1*log(seq(15,40)))

fit = nls(dyld ~ eval_scaled_point(tmx,deltat,coef0,coef1), data = fdata, start = list(coef0 = scoef0, coef1 = scoef1))
coef0 = fit$m$getPars()[["coef0"]]
coef1 = fit$m$getPars()[["coef1"]]

eval_scaled_point(26.6,6,coef0,coef1)

pdata$dyldnew = eval_scaled_point(pdata$tmx,pdata$deltat,28.8,2)

ggplot(pdata) + 
  geom_line(aes(x = deltat, y = dyldnew, color = as.factor(tmx)), size = 1) + ylim(-1,0.2)

pdata$dyldnew = eval_scaled_step_point(pdata$tmx,pdata$deltat)

ggplot(pdata) + 
  geom_line(aes(x = deltat, y = dyldnew, color = as.factor(tmx)), size = 1.5) + ylim(-1,0.2)


ggplot(pdata) + 
  geom_line(aes(x = deltat, y = dyld, color = as.factor(tmx)), size = 1) + ylim(-1,0.1)

k = 1
a = 0.39
m = 28.8
b = 2
x = seq(20,35,0.1)
scalesens = k - (k-a)/(1 + exp(-b*(x-m)))
plot(x,scalesens)

bhsenscoef0 = -0.009038992944479
bhsenscoef1 = 0.001231793519879

bhdata = data.frame(tmx = seq(20,40,0.1))
bhdata$baseedd = eval_nxdd(betasedd,bhdata$tmx)
bhdata$bh = bhsenscoef0 + bhsenscoef1*log(bhdata$baseedd)
ggplot(bhdata) + geom_line(aes(x = baseedd, y = bh))

ggplot(statdata,aes(x = tmx, y = eddsens)) +
  geom_point(aes(shape = modeltype, color = bootstrap.group),size = 4) + 
  geom_line(data = bhdata,aes(x = tmx,y=bh)) +
  theme_classic()

ggplot() + geom_line(data = bhdata,aes(x = tmx,y=bh))



#This is just reproducibility testing
# caldata = read.csv("rasters_sr_memo_global/Sacks_Zarc_fill_fill_120d/butler1/sr_deltaT2_2000.csv")
# caldata = filter(caldata,Crop == crop)
# caldata$tmx = caldata$tmaxbase
# caldata$tmp = caldata$tempbase
# # caldata$X.Temp = caldata$ludtemp1.00
# # caldata$deltay = caldata$noscaledyppDLU1.00
# caldata$X.Temp = caldata$ludtemp0.00
# caldata$deltay = caldata$noscaledyppDLU0.00
# 
# # caldata$deltagdd = caldata$dgddDLU1.00
# # caldata$deltaedd = caldata$deddDLU1.00
# 
# # plot(caldata$deltaedd,caldata$deddDLU0.00)
# # plot(caldata$deddDLU1.00,caldata$deltay)
# plot(caldata$tmx,caldata$eddsens)
# summary(caldata$eddsens)
