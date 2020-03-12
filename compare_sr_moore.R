# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("ggridges","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# Setup for SR
tdbase = "Sacks_Zarc_fill_fill_120d"
srinfolder = paste0("rasters_sr_memo_global/",tdbase,"/")
moinfolder = paste0("../AgroServYield_coefficients/estimated_global/Feb7/")

srscens = c("deltaT0","deltaT1","deltaT2")
srcrops = c("Maize","Soybeans","Cotton")

moscens = c("deltaTssp0","deltaT1degC","deltaT2degC")
mocrops = c("Maize","Rice","Wheat")

names(moscens) <- srscens

sryear = 2000
moyear = 2050

dlus = c(0.00,0.25,0.50,1.00)
dlusufs = sprintf("%.2f",dlus)

# File with the area weights for all crops
wgtfname = "GIS/areas.csv"

# File with equivalencies for conversion from the polygon grid to the points grid
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

# Zones
zonesfname = "GIS/COLROW30_K_G.shp"

# Read the SR results
srdata = data.frame()

for (iscen in srscens) {
  infname = paste0(srinfolder,"sr_",iscen,"_",sryear,".csv")
  indata = read.csv(infname)
  srdata = rbind(srdata,indata)
}

# Keep only relevant variables (SR) and gather by DLU
srdata = srdata[c("COLROW30","Crop","scen",paste0("dyppDLU",dlusufs))]
srdata = gather(srdata,dlu,dyld,paste0("dyppDLU",dlusufs),factor_key = T)
levels(srdata$dlu) <- regmatches(levels(srdata$dlu),regexpr("[0-9].[0-9][0-9]",levels(srdata$dlu)))

# Build the Moore crops with SR estimates
tempcot = filter(srdata,Crop == "Cotton")
names(tempcot)[names(tempcot) == "dyld"] <- "V1"
tempcot$Crop <- NULL
tempsoy = filter(srdata,Crop == "Soybeans")
names(tempsoy)[names(tempsoy) == "dyld"] <- "V2"
tempsoy$Crop <- NULL
temp = left_join(tempcot,tempsoy)
temp$dyld = (temp$V1 + temp$V2)/2.0
temp$V1 <- NULL
temp$V2 <- NULL
temp$Crop = "Rice"
srdata = rbind(srdata,temp)
temp$Crop = "Wheat"
srdata = rbind(srdata,temp)

srdata$method = "SR"

# Read the Moore results
modata = data.frame()

for (iscen in moscens) {
  infname = paste0(moinfolder,iscen,"/estimate_range_",iscen,"_",moyear,"_all.rds")
  indata = readRDS(infname)
  modata = rbind(modata,indata)
}

#Filter by year
modata = filter(modata,ScenYear == moyear)

# Keep only relevant variables (MO) and gather by DLU
names(modata)[names(modata) == "ID"] <- "COLROW30"
names(modata)
modata = modata[c("COLROW30","Crop","scen",paste0("dlu",dlusufs))]
modata = gather(modata,dlu,dyld,paste0("dlu",dlusufs),factor_key = T)
levels(modata$dlu) <- regmatches(levels(modata$dlu),regexpr("[0-9].[0-9][0-9]",levels(modata$dlu)))

# Build the SR crops with Moore estimates
tempric = filter(modata,Crop == "Rice")
names(tempric)[names(tempric) == "dyld"] <- "V1"
tempric$Crop <- NULL
tempwht = filter(modata,Crop == "Wheat")
names(tempwht)[names(tempwht) == "dyld"] <- "V2"
tempwht$Crop <- NULL
temp = left_join(tempric,tempwht)
temp$dyld = (temp$V1 + temp$V2)/2.0
temp$V1 <- NULL
temp$V2 <- NULL
temp$Crop = "Soybeans"
modata = rbind(modata,temp)
temp$Crop = "Cotton"
modata = rbind(modata,temp)

modata$method = "Moore"


# Bind both methods
alldata = rbind(srdata,modata)

# Read fractional area weights
wgtdata = read.csv(wgtfname)
wgtdata = gather(wgtdata,Crop,Area,Wheat:Cotton,factor_key = T)

# Read conversion names ID->COLROW30
convdata = read.csv(convfname)

wgtdata = left_join(wgtdata,convdata, by = c("ID" = "final.ID"))

alldata = left_join(alldata,wgtdata, by = c("COLROW30" = "final.COLROW30", "Crop" = "Crop"))

zonesshp = st_read(zonesfname)
zonesshp$zonename <- NA
zonesshp$zonecode <- NA
zonesshp$zonename[zonesshp$Class4 == "warm temperate"] <- "TEMP"
zonesshp$zonename[zonesshp$Class4 == "equatorial"] <- "EQUA"
zonesshp$zonename[zonesshp$Class4 == "arid"] <- "ARID"
zonesshp$zonename[is.na(zonesshp$zonename)] <- "OTHR"
zonesshp$zonecode[zonesshp$Class4 == "warm temperate"] <- 1
zonesshp$zonecode[zonesshp$Class4 == "equatorial"] <- 2
zonesshp$zonecode[zonesshp$Class4 == "arid"] <- 3
zonesshp$zonecode[is.na(zonesshp$zonecode)] <- 4
zonecodes = c(1,2,3,4)
names(zonecodes) = c("TEMP","EQUA","ARID","OTHR")

zonesdata = zonesshp[c("COLROW30","zonecode","zonename")]
zonesdata$geometry <-NULL


alldata = left_join(alldata,zonesdata, by = c("COLROW30"))
alldata$Area[is.na(alldata$Area)] <- 0.0

# Sync names of scenarios
for (i in 1:length(moscens)) {
  iscen = moscens[i]
  alldata$scen[alldata$scen == iscen] <- names(iscen)
}

# Convert some factor variablx`es that ended up becoming characters
alldata$COLROW30 <- as.factor(alldata$COLROW30)
alldata$ID <- as.factor(alldata$ID)
alldata$Country <- as.factor(alldata$Country)
alldata$method <- as.factor(alldata$method)
alldata$Crop <- as.factor(alldata$Crop)
alldata$zonename <- as.factor(alldata$zonename)
alldata$zonecode <- as.factor(alldata$zonecode)

alldata$scen = droplevels(alldata$scen)

# shpsrdata = left_join(zonesshp[c("geometry","COLROW30")],srdata, by = "COLROW30")
# plot(filter(shpsrdata, Crop == "Soybeans" & scen == "deltaT0" & Area >= 0.001)[c("geometry","Area")], pch = 20)
# write.csv(alldata,"alldata.csv", row.names = F)

# Clean house, this takes up a lot of RAM
rm(temp,tempcot,tempric,tempwht,tempsoy)
rm(modata,srdata)


###############################################################
repdata = alldata
repdata <- filter(alldata,Area >= 0.001)
repdata$reps = ceiling(repdata$Area/0.002)
repdata <- repdata[rep(row.names(repdata), repdata$reps), ]
# plot(ecdf(repdata$Area), xlim = c(0.0,0.01))

crop = "Soybeans"
# FIXME doesn't work passing crop to filter
cropdata <- alldata %>% filter(Crop == crop)
repcropdata <- repdata %>% filter(Crop == crop)


ggplot(filter(repcropdata, zonename != "OTHR"),
       aes(x = dyld, y = dlu)) + 
  facet_grid(scen + zonename ~ method) +
  geom_density_ridges(aes(fill = zonename, alpha = 0.8), show.legend = F) +
  theme_ridges() + 
  labs(x = "Yield change (%)",
       y = "Land use change (pp of whole cell)",
       title = crop
  )

ggplot(filter(repcropdata, Country == "Brazil" & zonename != "OTHR"),
       aes(x = dyld, y = dlu)) + 
  facet_grid(scen + zonename ~ method) +
  geom_density_ridges(aes(fill = zonename, alpha = 0.8), 
                      show.legend = F) +
  theme_ridges() + 
  labs(x = "Yield change (%)",
       y = "Land use change (pp of whole cell)",
       title = crop,
       subtitle = "Brazil"
  )



for (crop in c("Soybeans","Maize","Wheat","Cotton","Rice")) {
repcropdata <- repdata %>% filter(Crop == crop)
plt <- ggplot(filter(repcropdata, Country == "Brazil" & zonename != "OTHR"),
       aes(x = dyld, y = dlu)) + 
  facet_grid(scen + zonename ~ method) +
  geom_density_ridges(aes(fill = zonename, alpha = 0.8), 
                      show.legend = F) +
  theme_ridges() + 
  labs(x = "Yield change (%)",
       y = "Land use change (pp of whole cell)",
       title = crop,
       subtitle = "Brazil"
  )
print(plt)
}

#Testing
ggplot(filter(repcropdata, Country == "Brazil" & zonename != "OTHR" & !(scen == "deltaT0" & dlu == "0.00")),
       aes(x = dyld, y = dlu)) + 
  facet_grid(scen + zonename~ method) +
  geom_density_ridges(aes(fill = zonename, alpha = 0.8),
                      show.legend = F) +
  theme_ridges() + 
  labs(x = "Yield change (%)",
       y = "Land use change (pp of whole cell)",
       title = crop,
       subtitle = "Brazil"
  ) + 
  geom_vline(xintercept = -5)

# 
# theme_Publication <- function(base_size=10, base_family="Helvetica") {
#   
#   library(grid)
#   
#   library(ggthemes)
#   
#   (theme_foundation(base_size=base_size, base_family=base_family)
#     
#     + theme(plot.title = element_text(face = "bold",
#                                       
#                                       size = rel(1.2), hjust = 0.5),
#             
#             text = element_text(),
#             
#             panel.background = element_rect(colour = NA),
#             
#             plot.background = element_rect(colour = NA),
#             
#             panel.border = element_rect(colour = NA),
#             
#             axis.title = element_text(face = "bold",size = rel(1)),
#             
#             axis.title.y = element_text(angle=90,vjust =2),
#             
#             axis.title.x = element_text(vjust = -0.2),
#             
#             axis.text = element_text(), 
#             
#             axis.line = element_line(colour="black"),
#             
#             axis.ticks = element_line(),
#             
#             panel.grid.major = element_blank(),
#             
#             panel.grid.minor = element_blank(),
#             
#             legend.key = element_rect(colour = NA),
#             
#             legend.position = "right",
#             
#             legend.direction = "horizontal",
#             
#             legend.key.size= unit(0.4, "cm"),
#             
#             legend.margin = margin(t = c(3,3,3,3), unit='mm'),
#             
#             legend.title = element_text(face="italic"),
#             
#             plot.margin=margin(t = c(5,5,5,5), unit='mm'),
#             
#             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#             
#             strip.text = element_text(face="bold")
#             
#     ))
#   
#   
#   
# }