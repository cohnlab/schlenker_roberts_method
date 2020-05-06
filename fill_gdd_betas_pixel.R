library(tidyverse)
library(distances)

calname = "Sacks_ZARC_fill_fill_120d"
ddversionstring = "rgagcfsr_bound_1995_2005"
infolder = paste0("gdd_betas_pixel/",calname,"/",ddversionstring,"/")
outfolder = paste0("gdd_betas_pixel_fill/",calname,"/",ddversionstring,"/")

idvname = "COLROW30"
dvnames = c("beta1","beta2","beta3")
drvname = "beta1"

# Number of nearest neighbors to check for values to fill
# Ideally should be the number of NAs + 1, but setting to lower improves performance
knn = 10000

shpfname = "GIS/COLROW30.shp"

shp = st_read(shpfname) 

ddata = shp %>% st_set_geometry(NULL) %>% dplyr::select(-RealArea_m)

allids = ddata %>% dplyr::select(COLROW30)

dobj = distances(ddata,id_variable = "COLROW30")

infsufs = dir(infolder)

dir.create(outfolder, showWarnings = F, recursive = T)

for (infsuf in infsufs) {
  print(infsuf)
  infname = paste0(infolder,infsuf)
  outfname = paste0(outfolder,infsuf)
  
  indata = read.csv(infname)
  alldata = allids %>% left_join(indata)
  filldata = alldata
  
  idvcol = which(names(alldata) %in% idvname)
  dvcols = which(names(alldata) %in% dvnames)
  drvcol = which(names(alldata) %in% drvname)
  
  nainds = which(is.na(alldata[drvname]))
  naids = alldata[nainds,idvcol]
  
  # i = 100
  pb <- txtProgressBar(min = 0, max = length(nainds), style = 3)
  for (i in 1:length(nainds)) {
    # print(i)
    ind = nainds[i]
    nninds = nearest_neighbor_search(dobj,knn,query_indices = ind)
    nnind = nninds[min(which(!is.na(alldata[nninds,drvcol])))]
    filldata[ind,dvcols] <- alldata[nnind,dvcols]
    setTxtProgressBar(pb, i)
  }
  
  write.csv(filldata,file = outfname,row.names = F)
}