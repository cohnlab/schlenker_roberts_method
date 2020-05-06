###########################################################################
#
#   This reads a set of coefficent tables 
#   and generates estimates of yield impact for various LU change levels
###########################################################################

# --- ---- Packages # --- ---- # --- ----
Packages <- c("dplyr","tidyverse","data.table","sf",
              "ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer")
lapply(Packages, library, character.only = TRUE)

# --- ---- # --- ---- # --- ---- # --- ----


# --- ---- Folders # --- ----# --- ----
# infolder <- "C:/Users/u244034/Dropbox/AgroServ/Work/Agroserv_Yield/AgroServYield_coefficients/Result analysis/Inputs/"
# outfolder <- "C:/Users/u244034/Dropbox/AgroServ/Work/Agroserv_Yield/AgroServYield_coefficients/Result analysis/Outputs/"
infolder <- "deploy_tests/Coefficient Tables/"
outfolder <- "deploy_tests/"



MRfolder <- "MR approach/"
AGTfolder <- "AgroservTemp/"

# --- ----# --- ----# --- ----# --- ----


# --- ---- Files # --- ----# --- ----# --- ----

adpfname = paste0(infolder,MRfolder,"/AD.csv") # Adapt coef values
betfname = paste0(infolder,MRfolder,"/Betas.csv") #Moore's betas
tbsfname = paste0(infolder,MRfolder,"/Tbaseline.csv") # Tbaseline for each gridcell


# dttfname = paste0(infolder,AGTfolder,"/deltaTssp245.csv") # Tbaseline for each gridcell
dttfname = paste0(infolder,AGTfolder,"/deltaTrcp245.csv") # Tbaseline for each gridcell
agtfname.nl = paste0(infolder,AGTfolder,"/NL.csv") # AgroservT coefficient values
agtfname.l = paste0(infolder,AGTfolder,"/L.csv") # AgroservT coefficient values


# --- ---- Other parameters # --- ----# --- ----

crops = c("Maize","Rice","Wheat")

syst = "HI"
year = 2050

# The values of change in land use to evaluate. Fractional

dlus.l = c(0.00,0.25,0.5,0.75,1.0)
dlus.nl = c(0.00,0.25,0.5,0.75,1.0)

scenarios.ssp = c("deltaTrcp245")
scenarios.t = c("deltaT0","deltaT1","deltaT2","deltaT3","deltaT4","deltaT5")


# --- ----# OPEN GENERAL FILES  # --- ----# --- ----# --- ----# --- ----

adpdata <- read.csv(adpfname) %>% filter(MngSystem == syst) 
betdata <- read.csv(betfname)
tbsdata <- read.csv(tbsfname)
dttdata <- read.csv(dttfname)
agtdata.nl <- read.csv(agtfname.nl)
agtdata.l <- read.csv(agtfname.l)


# Long to wide
adpdata <- spread(adpdata,"PARAMETER","Value")
betdata <- spread(betdata,"PARAMETER","Value")
tbsdata <- spread(tbsdata,"PARAMETER","Value")
dttdata <- spread(dttdata,"PARAMETER","Value")
agtdata.nl <- spread(agtdata.nl,"PARAMETER","Value")
agtdata.l <- spread(agtdata.l,"PARAMETER","Value")

agtdata.nl$NL <- agtdata.nl$NL/2
agtdata.l$L <- agtdata.l$L/2

# Join them
alldata <- left_join(dttdata,tbsdata, by = c("ID","Crop")) %>% 
  left_join(betdata, by = c("Crop")) %>%
  left_join(agtdata.nl, by = c("ID")) %>%
  left_join(agtdata.l, by = c("ID"))
alldata$AD <- adpdata$AD
alldata$ID <- as.factor(alldata$ID)

alldata = filter(alldata,ScenYear %in% year)
alldata = filter(alldata,Crop %in% crops)
alldata = subset(alldata, select = -gccdtmax )

#########   write standardized Moore table #########

write.csv(alldata,paste0(outfolder,"Moore_table.csv"))

######### ######### ######### ######### ######### #########

tempdata <- subset(alldata, select = c(ID,Crop,scen,gccdtemp))
tempdata <- spread(tempdata,"scen","gccdtemp")

for (scen in scenarios.t){
  tempdata[,scen] <- as.integer(substr(scen,nchar(scen),nchar(scen)))}


newdata <- subset(alldata, select = -c(scen,ScenYear,gccdtemp) )


compute_agroservY <- function(col,dlu1,dlu2) {
  with(newdata,
       beta_1*(col+(L*dlu1) + (NL*dlu2)) + 
         beta_2*(col+(L*dlu1) + (NL*dlu2))^2 + 
         beta_3*Tbaseline*(col+(L*dlu1) + (NL*dlu2)) + 
         beta_4*Tbaseline*(col+(L*dlu1) + (NL*dlu2))^2 + 
         beta_5*AD*(col+(L*dlu1) + (NL*dlu2)))
}

compute_agroservDeltaT <- function(col,dlu1,dlu2) {
  with(newdata,col + (L*dlu1) + (NL*dlu2))
}



scenarios <- c("ssp245",scenarios.t)

  LUdy <- do.call(rbind,
            lapply(scenarios, 
                function(sc) do.call(rbind,
                  lapply(dlus.l, 
                    function(x) do.call(rbind,
                      lapply(dlus.nl, function(y) 
                          data.frame(tempdata[,1:2],
                               dyld=compute_agroservY(tempdata[,sc],x,y),
                               dtemp=compute_agroservDeltaT(tempdata[,sc],x,y),
                               dluNL=y,
                               dluL=x,
                               Scen=sc)))))))
  LUdy$dtmax = LUdy$dtemp
  LUdy = LUdy %>% 
    dplyr::select(ID,Scen,Crop,dluL,dluNL,dyld,dtemp,dtmax) %>% 
    dplyr::rename(scen = "Scen")
  
  saveRDS(LUdy,paste0(outfolder,"results_table_MR.rds"))
    
# for (ccrop in crops)
# {
#   final.table <- filter(LUdy,Crop %in% ccrop)
#   write.csv(LUdy,paste0(outfolder,ccrop,".csv"))
# }
  
