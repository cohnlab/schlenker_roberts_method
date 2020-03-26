require(tidyverse)

infolder = "../globiom_clim_extract/output_tables/"
outfolder = "scenario_tables/"

invarnames = c("gccdtemp","gccdtmax")
outvarnames = c("gccdtemp","gccdtmax")
names(outvarnames) <- invarnames

crops = c("Maize","Soybeans","Cotton")
scenarios = c("ssp245","ssp585")

# Create output folder
dir.create(outfolder,showWarnings = F, recursive = T)

alldata = data.frame()

# crop = crops[1]
# scenario = scenarios[1]
# invarname = invarnames[1]
for (crop in crops) {
  for (scenario in scenarios) {
    for (invarname in invarnames) {
      outvarname = outvarnames[[invarname]]
      
      infname = paste0(infolder,crop,".",scenario,".",invarname,".all.csv")
      indata = read.csv(infname)
      indata %>% 
        dplyr::select(-c("XCOORD","YCOORD","RealArea_m")) %>% 
        rename(ID = "COLROW30", Crop = "crop", scen = "rcp") %>%
        gather(key = "ScenYear",value = "Value",-ID,-Crop,-scen) -> indata
      indata$ScenYear = as.factor(as.numeric(substr(indata$ScenYear,2,5)))
      indata$PARAMETER = outvarname
      
      alldata = rbind(alldata,indata)
    }
  }
}
# Put zero where NA
alldata$Value[is.na(alldata$Value)] <- 0

for (scenario in scenarios) {
  outfname = paste0(outfolder,"deltaT",scenario,".csv")
  alldata %>% filter(scen == scenario) -> outdata
  write.csv(outdata, file = outfname, row.names = F)
}