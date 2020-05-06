Packages <- c("Rcpp","dplyr","tidyverse","gridExtra","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

# Function that applies the DD~T curves
source("eval_nxdd.R")

tdbase = "Sacks_ZARC_fill_fill_120d"

# ddversionstring = "agcfsr_unbound_usa"
ddversionstring = "agcfsr_bound_usa"
# ddversionstring = "merra"
# ddversionstring = "wgt_crop"
# ddversionstring = "wgt_crop_un1bound"

crops = c("Maize","Soybeans","Cotton")
# crops = c("Maize","Soybeans")

# Base sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Reference countries and what impacts should be at hirefcountry's Tmax at reference deltat. 
hirefcountry = "BRA"
lorefcountry = "USA"
refdeltat = 6.0
# refmult = 1.25
refmult = 1.50
refvals = c(-55.88, -48.82, -46.72)*refmult
names(refvals) <- crops
# Use computed instead of the reported refvals
icomputerefvals = TRUE

# XDD ~ Tx functions
betasfolder = paste0("gdd_betas/",tdbase,"/",ddversionstring,"/")

# Base data to evaluate functions on (tmx and tmp are Tbaselines)
bigbasedata = readRDS("tbaseline_countries.rds")

# Output file name
outfolder = paste0("tmaxsens/",ddversionstring,"/")
outfname = paste0(outfolder,"tmaxsens_linear_2_",refmult,".csv")

# SR values to compare to
srrefdata = read.csv("aux_figures/test_table_A5.csv") %>% filter(source == "SR Table A5")

# GLOBIOM grid, points shapefile and reference for rasterization
shpfname = "GIS/COLROW30.shp"
cntfname = 'GIS/ne_50m_admin_0_countries_lakes.shp'
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

shp = st_read(shpfname)
cntshp = st_read(cntfname)
st_crs(shp) <- st_crs(cntshp)
shp = st_join(shp,cntshp)[c("COLROW30","ADM0_A3")]

refrast = raster(reffname)

# Plot maps of impacts on the GLOBIOM grid
iglobplot = TRUE

# Base temperatures at the GLOBIOM grid
tempbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmp.csv")
tmaxbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmx.csv")

tempbase = read.csv(tempbasefname) %>% gather("Crop","tempbase",-ID)
tmaxbase = read.csv(tmaxbasefname) %>% gather("Crop","tmaxbase",-ID)
tempdata = left_join(tmaxbase,tempbase, by = c("ID","Crop"))

areadata = read.csv("GIS/areas.csv")
shp = left_join(shp,areadata, by = "COLROW30")
shp$COLROW30 = as.factor(shp$COLROW30)

dir.create(outfolder,showWarnings = F, recursive = T)

# caldata = data.frame(tmx = c(27,27,30), X.Temp = c(1,6,1), deltay = c(-6.4,-55.9,-1))
# caldata = read.csv("../moore_coefficients/includedstudies_revisions_tbaseline.csv")
# names(caldata)[names(caldata) == "Yield.change...."] <- "deltay"


# Functions to evaluate the method. Depend on a lot of variables being in .GlobalEnv
# Evaluate using a fixed scalesens. Set it to 1 to get the original value
eval_fixscale_point <- function(tmx,tmp,deltat,scalesens) {
  # tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
  # return(tmp)
}
# Evaluate using a stepwise tmaxbase~scalesens relationship
eval_scaled_step_point <- function(tmx,deltat,coeftab) {
  tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  scalesens = coeftab$inter[1] + coeftab$slp[1]*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(coeftab$minscale[1],scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}

# Evaluate using a stepwise tmaxbase~scalesens relationship, taking tmp as well
eval_scaled_step_point_tmp <- function(tmx,tmp,deltat,cropcoeftab) {
  # tmp = fittmx$coefficients[[1]] + tmx*fittmx$coefficients[[2]]
  dgdd = eval_nxdd(betasgdd,tmp+deltat) - eval_nxdd(betasgdd,tmp)
  dedd = eval_nxdd(betasedd,tmx+deltat) - eval_nxdd(betasedd,tmx)
  
  scalesens = cropcoeftab$inter[1] + cropcoeftab$slp[1]*tmx 
  scalesens = pmin(1,scalesens)
  scalesens = pmax(cropcoeftab$minscale[1],scalesens)
  dyld = dgdd*gddsens[[crop]] + dedd*eddsens[[crop]]*scalesens
  return(dyld)
}


coeftab = data.frame()
repdata = data.frame()
plots_range = list()
plots_fun = list()
croprasts = list()
# crop = "Maize"
for (crop in crops) {
  
  # Read the crop-specific coefficients for the T-GDD functions
  betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
  betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))
  
  # basedata per crop
  basedata = filter(bigbasedata,Crop == crop)
  
  # Get Tmean Tmax relationship
  # ggplot(basedata) + geom_point(aes(x = tmp, y = tmx))
  # summary(lm(tmx ~ tmp, basedata))
  # summary(lm(tmp ~ tmx, basedata))
  fittmx = lm(tmp ~ tmx, basedata)
  
  hireftmx = filter(basedata,country == hirefcountry)$tmx
  loreftmx = filter(basedata,country == lorefcountry)$tmx
  hireftmp = filter(basedata,country == hirefcountry)$tmp
  loreftmp = filter(basedata,country == lorefcountry)$tmp
  
  # Compute loref values to compare with reported data
  dumrefdata = data.frame(deltat = seq(1,10),
                          dyld = eval_fixscale_point(loreftmx,loreftmp,seq(1,10),1)*100.0,
                          Crop = crop, source = "Average Tbaselines")
  srrefdata = rbind(srrefdata,dumrefdata)
  
  # Use computed or reported refvals (either must be in percentage)
  if (icomputerefvals) {
    refval = eval_fixscale_point(loreftmx,loreftmp,refdeltat,1)*refmult*100.0
  } else {
    refval = refvals[[crop]]
  }
  
  
  # Find at which scalesens impacts at hireftmx are refval
  # Find the root of the function - refval
  rootobj = uniroot(function(x) {eval_fixscale_point(hireftmx,hireftmp,refdeltat,x)} - refval*0.01,
                    c(0,1))
  minscale = rootobj$root
  
  # Get the intercept and slope of the line
  fit = lm(scalesens ~ tmx,
           data = data.frame(tmx = c(loreftmx,hireftmx), scalesens = c(1.0,minscale)))
  
  cropcoeftab = data.frame(Crop = crop,
                           lotmx = loreftmx,
                           hitmx = hireftmx,
                           refval = refval,
                           inter = fit$coefficients[[1]],
                           slp = fit$coefficients[[2]],
                           minscale = minscale)
  coeftab = rbind(coeftab,cropcoeftab)
  
  
  # Synthetic data for plotting 
  pdata = data.frame(tmx = seq(15,36,3))
  pdata$tmp = fittmx$coefficients[[1]] + pdata$tmx*fittmx$coefficients[[2]]
  dumdata = pdata
  pdata = data.frame()
  for (dt in seq(0,10,0.5)) {
    dumdata$deltat = dt
    pdata = rbind(pdata,dumdata)
  }
  pdata$dgdd = eval_nxdd(betasgdd,pdata$tmp+pdata$deltat) - eval_nxdd(betasgdd,pdata$tmp)
  pdata$dedd = eval_nxdd(betasedd,pdata$tmx+pdata$deltat) - eval_nxdd(betasedd,pdata$tmx)
  
  pdata$dyld = eval_scaled_step_point(pdata$tmx,pdata$deltat,cropcoeftab)
  
  # Plot the estimated range of impacts
  plotrange = ggplot(pdata) + 
    geom_line(aes(x = deltat, y = dyld*100, color = as.factor(tmx)), size = 1.5) + 
    xlab(expression(paste(Delta,"T (°C)"))) +
    ylab(expression(paste(Delta,"T (pp)"))) +
    scale_colour_discrete(name="Tmax") +
    labs(title = crop) + ylim(-200,50)
  
  # Plot the actual scalesens function
  dumdata = data.frame(tmx = seq(15,40,0.1))
  dumdata$scalesens = cropcoeftab$inter[1] + cropcoeftab$slp[1]*dumdata$tmx
  dumdata$scalesens = pmin(1,dumdata$scalesens)
  dumdata$scalesens = pmax(cropcoeftab$minscale,dumdata$scalesens)
  plotfun = ggplot(dumdata) +
    geom_line(aes(x = tmx, y = scalesens), size = 1.5) +
    xlab("Tmax") +
    ylab("EDD sensitivity scale factor") +
    ylim(0.2,1.0) + 
    annotate(geom="text", x = 35, y = minscale+0.01, label = sprintf("%.2f",minscale), 
             size = 5, hjust = 0, vjust = 0) +
    annotate(geom = "text", x = 15, y = 0.6, size = 5, hjust = 0, vjust = 0, 
             label = paste0("US Tmax: ",sprintf("%.1f",loreftmx))) +
    annotate(geom = "text", x = 15, y = 0.55, size = 5, hjust = 0, vjust = 0, 
             label = paste0("US reference impact: ",sprintf("%.1f",refval/refmult),"%")) +
    annotate(geom = "text", x = 15, y = 0.4, size = 5, hjust = 0, vjust = 0, 
             label = paste0("BR Tmax: ",sprintf("%.1f",hireftmx))) +
    annotate(geom = "text", x = 15, y = 0.35, size = 5, hjust = 0, vjust = 0, 
             label = paste0("BR reference impact: ",sprintf("%.1f",refval),"%"))
    
  plots_range[[length(plots_range)+1]] = plotrange
  plots_fun[[length(plots_fun)+1]] = plotfun
  
  if (iglobplot) {
    cropshp = left_join(shp,filter(tempdata,Crop == crop),by = c("COLROW30" = "ID"))
    deltats = seq(1,6)
    scennames = paste0("deltaT",sprintf("%1i",deltats))
    for (i in 1:length(deltats)) {
      cropshp[scennames[i]] = eval_scaled_step_point_tmp(cropshp$tmaxbase,cropshp$tempbase,deltats[i],cropcoeftab)*100.0
    }
    valid = scennames
    croprast = rasterize(cropshp, refrast, field = valid, fun = mean, background = NA_real_,
                         by = NULL)
    croprasts[[crop]] <- croprast
    
    breaks = seq(-100,100,20)
    pal = brewer.pal(n = length(breaks), name = "RdBu")
    tit = crop
    
    plotobj <- tm_shape(croprast, bbox = c(-130,160,-50,65)) + 
      tm_raster(palette = pal, breaks = breaks,
                title = expression(paste(Delta,"Yield (%)"))) + 
      tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
                legend.outside = T) + 
      tm_layout(panel.show = T, panel.labels = scennames,
                panel.label.size = 1.2,
                main.title.position = "center", main.title = tit, main.title.size = 1.0) +
      tm_facets(nrow = 3, ncol = 2)
    print(plotobj)
    # ggsave(paste0(outfolder,"maps_",crop,".png"), plot = plotobj)
    
    uscropdata = cropshp %>% filter(ADM0_A3 == "USA") 
    
    for (i in 1:length(deltats)) {
      v = sum(uscropdata[[scennames[i]]]*uscropdata[[crop]],na.rm=T)/sum(uscropdata[[crop]],na.rm=T)
      vdata = data.frame(deltat = deltats[i], Crop = crop, dyld = v, source = "Impact average")
      repdata = rbind(vdata,repdata)
    }
  }
}

# Plot of composite C3 crops
if (iglobplot) {
  breaks = seq(-100,100,20)
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  tit = "Composite C3 (average of Soybeans and Cotton)"
  
  croprast = (croprasts$Soybeans + croprasts$Cotton )/2.0
  
  plotobj <- tm_shape(croprast, bbox = c(-130,160,-50,65)) + 
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = scennames,
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = 3, ncol = 2)
  print(plotobj)
  
  
  
  # write.csv(coeftab, file = outfname, row.names = F)
}

allrepdata = rbind(srrefdata,repdata)
allrepdata %>% 
  filter(source !=  "Average Tbaselines") %>%
  ggplot() + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source), size = 2) +
  xlim(1,6) +
  ylim(-100,10) +
  ggtitle("U.S. Area weighted impact averages", subtitle = ddversionstring)

allrepdata %>% spread(source,dyld) %>% 
  mutate(error = 100.0*(`Impact average` - `SR Table A5`)/`SR Table A5`) %>%
  ggplot() + geom_line(aes(x=deltat,y=error,color=Crop),size=2) +
  xlim(0,6) +
  ggtitle("Error in the U.S.", subtitle = ddversionstring)

allrepdata %>% 
  filter(source %in%  c("Average Tbaselines","SR Table A5")) %>%
  ggplot() + 
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source), size = 2) +
  xlim(1,6) +
  ylim(-60,10) +
  ggtitle("U.S.average temperature impacts", subtitle = ddversionstring)
  
# plots_all = list()
# for (i in seq(1,(2*length(plots_range)),2)) {
#   plots_all[i] = plot
# }

do.call(grid.arrange,c(plots_range,plots_fun,c(nrow = 2)))

ggplot(srrefdata) +
  geom_line(aes(x = deltat, y = dyld, color = Crop, linetype = source), size = 1.5)

# grid.arrange(plotrange,plotfun,nrow = 1)
