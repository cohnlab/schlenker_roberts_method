# Fits a DD~deltat model to each point

library(tidyverse)
library(raster)
library(sf)
library(broom)

calname = "Sacks_ZARC_fill_fill_120d"
# ddversionstring = "bcagcfsr"
ddversionstring = "rgagcfsr"

infolder  = paste0("rgagcfsr_dfs_partial/",calname,"/rgagcfsr/")
outfolder = paste0("gdd_betas_pixel_partial/",calname,"/",ddversionstring,"/")
# infolder  = "xavier_computed/"
# outfolder  = "xavier_dfs/"

dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

# crops = c("Soybeans")
# years = 2002:2002

crops = c("Maize","Cotton")
# crops = c("Maize","Soybeans","Cotton")
years = 2002:2008
# years = 2002:2002
deltats = 0:5 # Must get zero

# Points shapefile
shpfname = 'GIS/COLROW30.shp'
cntfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'

shp = st_read(shpfname)
cntshp = st_read(cntfname)
st_crs(shp) <- st_crs(cntshp)
cntdata = st_join(shp,cntshp)[c("COLROW30","ADM0_A3")] %>% st_set_geometry(NULL)

# crop = "Soybeans"
for (crop in crops) {
  print(crop)
  
  infname = paste0(infolder,crop,".climdata.deltat.rds")
  
  # Filter out points that have any NAs. We'll have to fill the coefficients later
  indata = readRDS(infname) %>% 
    filter(year %in% years) %>% 
    group_by(COLROW30) %>% 
    filter(all(!is.na(tempmean))) %>% 
    ungroup %>% as.data.frame()
  indata = left_join(indata,cntdata,by = c("COLROW30"))
  
  
  data = indata %>% 
    dplyr::select(c("COLROW30","YCOORD","XCOORD","ADM0_A3","year",
                    "deltat","agdd","bgdd","cgdd","fgdd"))
  
  # data = filter(data,ADM0_A3 %in% c("BRA","USA","CHI","RUS")) %>% 
  # group_by(ADM0_A3) %>% sample_n(2) 
  # data = filter(data,COLROW30 %in% c("252 - 212", "253 - 212"))
  
  # data = filter(data,ADM0_A3 %in% c("BRA","USA","CHI")) %>% droplevels
  # selpoints = sample(levels(data$COLROW30),12)
  # data = filter(data,COLROW30 %in% selpoints) %>% droplevels
  
  vnames = c("agdd","bgdd","cgdd","fgdd")
  # Calculate deltas
  for (vname in vnames) {
    tmp <- left_join(data,filter(data, deltat == 0)[c("COLROW30","year",vname)], by = c("COLROW30","year"))
    data[paste0("d",vname)] = tmp[paste0(vname,".x")] - tmp[paste0(vname,".y")]
  }
  
  for (vname in vnames) {
    # Fitting and unnesting with the coefficients
    print(paste0("Fitting ",vname,"..."))
    fm = as.formula(paste0("d",vname," ~ deltat + I(deltat^2) + I(deltat^3) - 1"))
    data %>% nest(-COLROW30) %>%
      mutate(fit = map(data, ~ lm(fm , data = .)),
             results = map(fit,augment),
             rsquared = map_dbl(map(fit,glance),"r.squared"),
             maxres = map_dbl(map(fit,residuals),max),
             coefs = map(fit,"coefficients")) %>% 
      unnest_wider(coefs) %>%
      dplyr::rename(beta1 = "deltat", beta2 = "I(deltat^2)", beta3 = "I(deltat^3)") -> vfits
    
    # Write full fit tibble as RDS and coefficient tables as CSV
    saveRDS(vfits, file = paste0(outfolder,crop,".fits.",vname,".rds"))
    vbetas = vfits %>% 
      dplyr::select(COLROW30,beta1,beta2,beta3) %>% 
      as.data.frame()
    write.csv(vbetas, file = paste0(outfolder,crop,".betas.",vname,".csv"), row.names = F)
    
    # Remove tibble to save memory
    rm(vfits)
  }
  
} #crops


# data %>% nest(-COLROW30) %>% 
#   mutate(fit = map(data, ~ lm(dEDD ~ deltat + I(deltat^2) + I(deltat^3), data = .)),
#          results = map(fit,augment),
#          coefs = map(fit,coefficients)) %>% unnest_wider(coefs)
# 
# data %>% nest(-COLROW30) %>% 
#   mutate(fit = map(data, ~ lm(dEDD ~ deltat + I(deltat^2) + I(deltat^3), data = .)),
#          results = map(fit,augment),
#          rsquared = map_dbl(map(fit,glance),"r.squared"),
#          coefs = map(fit,"coefficients")) 
# 
#   
# data %>% nest(-COLROW30) %>%
#   mutate(fit = map(data, ~ lm(fm, data = .)),
#          results = map(fit,augment)) %>% unnest(c(results)) %>%
#   ggplot() + geom_point(aes(x = deltat, y = dagdd)) +
#   geom_line(aes(x = deltat, y = .fitted)) +
#   facet_wrap(~COLROW30)
# 
# data %>% nest(-COLROW30) %>% 
#   mutate(fit = map(data, ~ lm(dGDD ~ deltat + I(deltat^2) + I(deltat^3), data = .)),
#          results = map(fit,augment)) %>% unnest(c(results)) %>%
#   ggplot() + geom_point(aes(x = deltat, y = dGDD)) + 
#   geom_line(aes(x = deltat, y = .fitted)) +
#   facet_wrap(~COLROW30)
# 
# data %>% nest(-COLROW30) %>%
#   mutate(fit = map(data, ~ lm(dEDD ~ deltat + I(deltat^2) + I(deltat^3) - 1 , data = .)),
#          results = map(fit,augment),
#          rsquared = map_dbl(map(fit,glance),"r.squared"),
#          maxres = map_dbl(map(fit,residuals),max),
#          coefs = map(fit,"coefficients")) %>% 
#   unnest_wider(coefs) %>%
#   dplyr::rename(beta1 = "deltat", beta2 = "I(deltat^2)", beta3 = "I(deltat^3)") -> fitsEDD
# 
# 
# data %>% nest(-COLROW30) %>%
#   mutate(fit = map(data, ~ lm(dEDD ~ deltat + I(deltat^2) + I(deltat^3) - 1 , data = .)),
#          results = map(fit,augment),
#          rsquared = map_dbl(map(fit,glance),"r.squared"),
#          coefs = map(fit,"coefficients")) -> big
# 
# ggplot(data) + geom_point(aes(x = deltat, y = dEDD, color = as.factor(year))) +
#   facet_wrap(~COLROW30)
# ggplot(data) + geom_point(aes(x = deltat, y = dGDD, color = as.factor(year))) +
#   facet_wrap(~COLROW30)
# 
# data %>% nest(-COLROW30) %>% mutate(fit = map(data, ~ lm(dEDD ~ deltat, data = .)),
#                                     results = map(fit,glance)) %>% unnest(results)
# 
# data %>% nest(-COLROW30) %>% mutate(fit = map(data, ~ lm(dEDD ~ deltat, data = .)),
#                                     results = map(fit,"coefficients")) %>% unnest(results)
# 
# data %>% nest(-COLROW30) %>% 
#   mutate(fit = map(data, ~ lm(dEDD ~ deltat + I(deltat^2) + I(deltat^3), data = .)),
#          results = map(fit,augment)) %>% mutate(fitted = map(results,.$.fitted))
# 
# rshp = left_join(shp,big[c("COLROW30","rsquared")])
# plot(rshp[c("rsquared","geometry")])

# 
# vname = "EDD"
# data = bind_rows(
#   by(data,droplevels(data$COLROW30),
#      function(d) {d[paste0("d",vname)] <- d[vname] - filter(d, deltat == 0)[[vname]];
#      return(d)}, 
#      simplify = F))
# 
# crop ="Soybeans"
# for (crop in crops) {
# print(crop)
#   readRDS(paste0(betasfolder,crop,".fitsGDD.rds")) %>% 
#   dplyr::select(COLROW30,beta1,beta2,beta3) %>% 
#   as.data.frame() %>%
#   write.csv(file = paste0(betasfolder,crop,".betasGDD.csv"), row.names = F)
# readRDS(paste0(betasfolder,crop,".fitsEDD.rds")) %>% 
#   dplyr::select(COLROW30,beta1,beta2,beta3) %>% 
#   as.data.frame() %>%
#   write.csv(file = paste0(betasfolder,crop,".betasEDD.csv"), row.names = F)
# }



