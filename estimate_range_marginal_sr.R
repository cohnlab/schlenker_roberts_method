# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer","classInt")
lapply(Packages, library, character.only = TRUE)

tdbase = "Sacks_Zarc_fill_fill_120d"

tempbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmp.csv")
tmaxbasefname = paste0("../AgroServYield_coefficients/Inputs_coef/",tdbase,"/Tbaseline_tmx.csv")

betasfolder = paste0("gdd_betas/",tdbase,"/")

# Points shapefile, reference raster and region of interest shapefile
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

# Mask, value is set depending on whether its global below
mskfname = "cropmap/cropland2000_grid_cru.nc"
imask = TRUE

iglobal = TRUE
if (iglobal) {
  regfname = "GIS/bbox_world.shp"
  globalbbox = c(-130,160,-50,65)
  # globalbbox = extent(c(230-360,300-360,25,50))
  
  msklim = 0.001
} else {
  # regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"
  regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"
  msklim = 0.000
}

brmskfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"

# Should we also plot CDFs?
iplotcdf = TRUE 
if (iplotcdf) {
  zonesfname = "GIS/COLROW30_K_G.shp"
}




# Path for output plots
if (iglobal) {
  outfolder = paste0("plots_sr_memo_global/",tdbase,"/")
  wrtfolder = paste0("rasters_sr_memo_global/",tdbase,"/")
} else {
  outfolder = paste0("plots_sr_memo/",tdbase,"/")
}

# crops = c("Soybeans")
# cropstrings = c("soybean")
crops = c("Maize","Soybeans","Cotton")
cropstrings = c("maize","soybean","cotton")
names(cropstrings) <- crops

# If true, use nGDD and nEDD scaling. Else apply coefficients directly as EDD and GDD
irescale = FALSE
if (irescale) {
# Number of days in the growing season used in SR20090
# ndays = c(180,180,210)
ndays = c(120,120,120)
names(ndays) <- crops
}

# Scale sensitivities based on VPD~EDD:Tmean relationships
ivpdscale = TRUE
if (ivpdscale) {
vpdscalefolder = paste0("vpd_scaling/",tdbase,"/")
vpdreftemp = 23 # Upper bound temp category of the reference region, should be the US here
}
 
# Cap for logY impacts, both positive and negative
icap = FALSE
cap = 0.5

# Sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Effect of 1pp of LU change in maximum and mean temperatures
izoneagt = TRUE
if (izoneagt) {
  agtfname = "../AgroServYield_coefficients/Resulting_tables/Feb7/agroservT.csv"
} else {
  luefftemp = 1.57*(0.5/0.7)
  luefftmax = 3.15*(0.5/0.7)
}

# FIXME: Fixed backgound warming.Substitute for a RCP file
bgtemp = 2.0
bgtmax = 2.0
basetit = paste0("BG + ",round(bgtemp),"\u00B0C")


# FIXME do for a single year for now
year = 2000

# Get a scenario string from the filename
# scen = tools::file_path_sans_ext(basename(dttfname))
scen = paste0("deltaT",bgtemp)

# Output plots
outfpref = paste0(outfolder,"/estimate_range_",scen,"_",year)


# Output plots
# outfpref = paste0(outfolder,"/estimate_range_",scen,"_",year)

# The values of change in land use to evaluate. Fractional
# dlus = c(0.00,1.00)
# dlus = c(1.00)
# dlus = c(0.0,0.3)
# dlus = c(0.0,0.3,1.00)
dlus = c(0.0,0.25,0.5,1.0)
# dlus = c(0.0)
# dlus = c(.00,0.05,0.10,0.20,0.3,0.5,1.0)
# dlus = c(0.00,0.05,0.10,0.20,0.3,1.0)

# Grid for arranging the plots
if (iglobal) {
  # gnrow = 2
  # gncol = ceiling(length(dlus)/gnrow)
  gnrow = length(dlus)
  gncol = 1
  # gncol = ceiling(length(dlus)/gnrow)
}else{
  gnrow = 1
  gncol = length(dlus)
}

# Maximum number of bins in plot. More than that, it increases bin size
maxbreaks = 30

# Override quantile based breaks anyway
ioverbreaks = TRUE
# overbreaks = seq(-30,30,5)
overbreaks = seq(-50,50,5)
# overbreaks = seq(-100,100,10)

# Create output folder
dir.create(outfolder, showWarnings = F)
dir.create(wrtfolder, showWarnings = F)

# Read baseline temp data
tempbase <- read.csv(tempbasefname) #%>% dplyr::select(-X)
tmaxbase = read.csv(tmaxbasefname) #%>% dplyr::select(-X)

# Make it long in respect to crops and join them
tempbase = gather(tempbase,"Crop","tempbase",-ID)
tmaxbase = gather(tmaxbase,"Crop","tmaxbase",-ID)
alldata <- left_join(tempbase,tmaxbase,by=c("ID","Crop"))

#FIXME: Just make a placeholder for GCC changes in temp and tmax
alldata$gccdtemp = bgtemp
alldata$gccdtmax = bgtmax

# Reag AGT coefficient for Tmean, and build ones for Tmax
if (izoneagt) {
  agtdata = read.csv(agtfname)
  agtdata <- spread(agtdata,"PARAMETER","Value")
  names(agtdata)[names(agtdata) == "agroservT"] <- "AGTtmp"
  agtdata$AGTtmx = agtdata$AGTtmp*2.0
  
  # Join them to alldata
  alldata <- left_join(alldata,agtdata,by=c("ID"))
}

# Using Rcpp here to avoid figuring out a faster lookup in R
# Evaluete a nEDD or nGDD function at a given T(max,mean)
cppFunction('NumericVector cpp_eval_nxdd(NumericMatrix x, NumericVector tvec) {
  int nrow = x.nrow(), ncol = x.ncol();
  int vsize = tvec.size();
  int p = 0;
  double t;
  double val;
  NumericVector vout(vsize);
  for (int p = 0; p < vsize; ++p) {
  t = tvec(p);
  for (int i = 0; i < nrow; ++i) {
  if (t >= x(i,0) && t<= x(i,1)) {
    vout(p) = x(i,2) + x(i,3)*t;
    break;
  }
  }
  }
  return vout;
}')

# Wrapper for cpp_eval_nxdd. Forces a single-column data.frame into a vector
eval_nxdd <- function(betas,tvec) {
  if (is.data.frame(tvec)) {
    tvec = tvec[,]
  }
  return(cpp_eval_nxdd(as.matrix(betas),as.vector(tvec)))
}



# Read the shapefile and the reference raster
shp = st_read(shpfname)
ref = raster(reffname)

# Read the relevant overlay shapefile, and set the bounding box to it if not global
reg = st_read(regfname)
if (iglobal) {
  bbox = extent(globalbbox)
} else {
  bbox = extent(reg)
}



# Read the mask file
if (imask) {
  msk = raster(mskfname)
  msk[msk<=msklim] <- NA
}

brshp = st_read(brmskfname)
brshp$msk = 1
brmsk = rasterize(brshp, ref, field = "msk", fun = mean, background = NA_real_,
                  by = NULL)

if (iglobal & imask) {
  msk[brmsk == 1] <- 1
}

if (iplotcdf) {
  zonesshp = st_read(zonesfname)
  zonesshp$zonename <- NA
  zonesshp$zonecode <- NA
  zonesshp$zonename[zonesshp$Class4 == "warm temperate"] <- "TEMP"
  zonesshp$zonename[zonesshp$Class4 == "equatorial"] <- "EQUA"
  zonesshp$zonename[zonesshp$Class4 == "arid"] <- "ARID"
  zonesshp$zonecode[zonesshp$Class4 == "warm temperate"] <- 1
  zonesshp$zonecode[zonesshp$Class4 == "equatorial"] <- 2
  zonesshp$zonecode[zonesshp$Class4 == "arid"] <- 3
  zonecodes = c(1,2,3)
  names(zonecodes) = c("TEMP","EQUA","ARID")
  
  zones = rasterize(zonesshp, ref, field = "zonecode", fun = mean, background = NA_real_,
                    by = NULL)
  
}



# Join all tables in a single file to be saved
outdata = data.frame()

# FIXME Loop this
# crop = "Soybeans"
for (crop in crops) {
  
  # Read the crop-specific coefficients for the T-GDD functions
  betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
  betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))
  
  # Filter baselines for the crop of interest
  cropdata <-alldata %>% filter(Crop ==  crop)
  
  # FIXME Filter the year, we should be able to loop this too
  # ydata <- cropdata %>% filter(ScenYear == year)
  ydata = cropdata
  
  # VPD scaling
  if (ivpdscale) {
    vpddata = read.csv(paste0(vpdscalefolder,crop,".vpdscale.csv"))
    cuts = vpddata$cuts
    vpddata$vpdscale = vpddata$coefs/vpddata$coefs[vpddata$cuts == vpdreftemp]
    ydata$tempcat = cut(ydata$tempbase,cuts)
    inds = findInterval(ydata$tempbase,cuts, left.open = T)
    inds[inds == 0] <- 1
    ydata$tempupper = vpddata$cuts[-1][inds]
    ydata$vpdscale = vpddata$vpdscale[-1][inds]
  }
  
  # Here we actually apply the equations, setting variable names prepended with the different LU levels
  ludtempvnames = paste0("ludtemp",sprintf("%.2f",dlus))
  ludtmaxvnames = paste0("ludtmax",sprintf("%.2f",dlus))
  ludgddvnames = paste0("dgddDLU",sprintf("%.2f",dlus))
  ludeddvnames = paste0("deddDLU",sprintf("%.2f",dlus))
  ludlogyvnames = paste0("dlogyDLU",sprintf("%.2f",dlus))
  lunoscaledlogyvnames = paste0("noscaledlogyDLU",sprintf("%.2f",dlus))
  ludyppvnames = paste0("dyppDLU",sprintf("%.2f",dlus))
  lunoscaledyppvnames = paste0("noscaledyppDLU",sprintf("%.2f",dlus))
  for (i in 1:length(dlus)) {
    # Changes in temp, tmax
    if (izoneagt) {
      ydata[ludtempvnames[i]] = ydata$AGTtmp*dlus[i]+ydata["gccdtemp"]
      ydata[ludtmaxvnames[i]] = ydata$AGTtmx*dlus[i]+ydata["gccdtmax"]
    } else {
      ydata[ludtempvnames[i]] = luefftemp*dlus[i]+ydata["gccdtemp"]
      ydata[ludtmaxvnames[i]] = luefftmax*dlus[i]+ydata["gccdtmax"]
    }
    
    # Changes in GDD, EDD, multiplying by ndays if necessary
    if (irescale) {
      ydata[ludgddvnames[i]] = (eval_nxdd(betasgdd,(ydata["tempbase"]+ydata[ludtempvnames[i]]) ) -
                                  eval_nxdd(betasgdd,ydata["tempbase"]))*ndays[crop]
      ydata[ludeddvnames[i]] = (eval_nxdd(betasedd,(ydata["tmaxbase"]+ydata[ludtmaxvnames[i]]) ) -
                                  eval_nxdd(betasedd,ydata["tmaxbase"]))*ndays[crop]
    } else {
      ydata[ludgddvnames[i]] = (eval_nxdd(betasgdd,(ydata["tempbase"]+ydata[ludtempvnames[i]]) ) -
                                  eval_nxdd(betasgdd,ydata["tempbase"]))
      ydata[ludeddvnames[i]] = (eval_nxdd(betasedd,(ydata["tmaxbase"]+ydata[ludtmaxvnames[i]]) ) -
                                  eval_nxdd(betasedd,ydata["tmaxbase"]))
    }                                  
    
    # Apply the GDD and EDD sensitivities

    
    # Apply VPD scaling only on EDD, not on GDD
    if (ivpdscale) {
      ydata[lunoscaledlogyvnames[i]] = gddsens[crop]*ydata[ludgddvnames[i]] +
        eddsens[crop]*ydata[ludeddvnames[i]]
      ydata[ludlogyvnames[i]] = gddsens[crop]*ydata[ludgddvnames[i]] +
        eddsens[crop]*ydata[ludeddvnames[i]]*ydata$vpdscale
    } else {
      ydata[ludlogyvnames[i]] = gddsens[crop]*ydata[ludgddvnames[i]] +
        eddsens[crop]*ydata[ludeddvnames[i]]
    }
    
    #FIXME: Cap impacts at cap
    if (icap) {
      ydata[[ludlogyvnames[i]]][ydata[[ludlogyvnames[i]]] > cap] = cap
      ydata[[ludlogyvnames[i]]][ydata[[ludlogyvnames[i]]] < -cap] = -cap
    }
    
    # Get a percent value
    ydata[ludyppvnames[i]] = ydata[ludlogyvnames[i]]*100.0
    if (ivpdscale) {
      ydata[lunoscaledyppvnames[i]] = ydata[lunoscaledlogyvnames[i]]*100.0
    }
  }
  
  # Join with the shapefile
  yshp = left_join(shp,ydata, by = c("COLROW30" = "ID"))
  
  # Save to output dataframe
  dumdata = yshp
  dumdata$geometry <- NULL
  dumdata$scen = scen
  dumdata$crop = crop
  outdata = rbind(outdata,dumdata)
  
  # yrast <- rasterize(yshp, ref, field = ludyppvnames, fun = "last", background = NA_real_,
                     # by = NULL)
  
  # Valid names to rasterize 
  valid = c("tempbase","tmaxbase",ludtempvnames,ludtmaxvnames,ludgddvnames ,ludeddvnames, ludlogyvnames, ludyppvnames, lunoscaledlogyvnames,lunoscaledyppvnames) 
  
  yrast <- rasterize(yshp, ref, field = valid, fun = mean, background = NA_real_,
                     by = NULL)
  
  writeRaster(yrast,paste0(wrtfolder,"/sr_",scen,"_",year,crop), overwrite = T)
  
  if (imask) {
    yrast = mask(yrast,msk)
  }
  
  # levs = c("dlu0.10","dlu0.20","dlu1.00")
  # levs = c("dlu0.05","dlu0.10","dlu0.20","dlu1.00")
  # subrast = subset(yrast,levs)
  subrast = subset(yrast,ludyppvnames)
  # subrast = subset(subrast,names(subrast)[-1])
  levs = names(subrast)
  
  qmin = round(quantile(as.array(subrast),0.05, na.rm=T))
  # breaks = seq(-5,5,0.5)
  breaks = seq(qmin,-1*qmin,2)
  if (length(breaks) > maxbreaks){
    breaks = seq(qmin,-1*qmin,5)
  }
  if (ioverbreaks) {
    breaks = overbreaks
  }
  
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  
  if (irescale) {
    tit = paste0("dYld | ndays = ",ndays[crop]," | ",crop, " | ",scen, " | ",tdbase," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  } else {
    tit = paste0("dYld | ",crop, " | ",scen, " | ",tdbase," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  }
  
  plotobj <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # Repeat without VPD scaling if that's enabled
  if (ivpdscale) {
    subrast = subset(yrast,lunoscaledyppvnames)
    
    if (irescale) {
      tit = paste0("Unscaled dYld", " | ",crop, " | ",scen," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
    } else {
      tit = paste0("Unscaled dYld | ", crop, " | ",scen," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
    }
    
    plotnoscale <- tm_shape(subrast, bbox = bbox) + 
      tm_raster(palette = pal, breaks = breaks,
                title = expression(paste(Delta,"Yield (%)"))) + 
      tm_shape(reg) + tm_borders() +
      tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
                legend.outside = T) + 
      tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
                panel.label.size = 1.2,
                main.title.position = "center", main.title = tit, main.title.size = 1.0) +
      tm_facets(nrow = gnrow, ncol = gncol)
  }
  
  # Repeating for the marginal effects
  subrast = subset(yrast,ludyppvnames)
  subrast = subrast - subrast$dyppDLU0.00
  names(subrast) <- ludyppvnames
  subrast = subset(subrast,names(subrast)[-1])
  levs = names(subrast)
  
  qmin = round(quantile(as.array(subrast),0.05, na.rm=T))
  # breaks = seq(-5,5,0.5)
  breaks = seq(qmin,-1*qmin,2)
  if (length(breaks) > maxbreaks){
    breaks = seq(qmin,-1*qmin,5)
  }
  if (ioverbreaks) {
    breaks = overbreaks
  }
  pal = brewer.pal(n = length(breaks), name = "RdBu")
  
  if (irescale) {
    tit = paste0("dYld (Marginal) | ndays = ",ndays[crop]," | ",crop, " | ",scen, " | ",tdbase," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  } else {
    tit = paste0("dYld (Marginal) | ",crop, " | ",scen, " | ",tdbase," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  }
  
  plotobjmarg <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = F, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # Other plots
  if (irescale) {
    titsuf = paste0("ndays = ",ndays[crop]," | ",crop, " | ",scen, " | ",tdbase)
  }else {
    titsuf = paste0(crop, " | ",scen, " | ",tdbase)
  }
  
  # dGDD
  subrast = subset(yrast,ludgddvnames)
  tit = paste0("dGDD | ",titsuf)
  breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
  plotdgdd <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = "YlOrRd", breaks = breaks,
              title = expression(paste(Delta,"GDD"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # dTmean
  subrast = subset(yrast,ludtempvnames)
  tit = paste0("dTmean | ",titsuf)
  breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
  plotdtmp <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = "YlOrRd", breaks = breaks,
              title = expression(paste(Delta,"Tmean"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  #dEDD
  subrast = subset(yrast,ludeddvnames)
  tit = paste0("dEDD | ",titsuf)
  breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
  plotdedd <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = "YlOrRd", breaks = breaks,
      title = expression(paste(Delta,"EDD"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # dTmax
  subrast = subset(yrast,ludtmaxvnames)
  tit = paste0("dTmax | ",titsuf)
  breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
  plotdtmx <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = "YlOrRd", breaks = breaks,
              title = expression(paste(Delta,"Tmax"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # Anomalies Scaled -  Unscaled
  subrast = subset(yrast,ludyppvnames) - subset(yrast,lunoscaledyppvnames)
  tit = paste0("Scaled - Unscaled")
  breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
  plotanom <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = "RdBu", midpoint = 0, breaks = breaks,
              title = expression(paste("Change in ",Delta,"Yield"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = T, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # CDFs
  if (iplotcdf) {
    png(filename = paste0(outfpref,crop,"_CDF_.png"), width = 600, height = 600, unit = "px", pointsize = 16)
    # Force a square set of plots here
    cdfnrow = ceiling((gnrow*gncol)/2)
    cdfncol = ceiling((gnrow*gncol)/2)
    par(mfrow = c(cdfnrow,cdfncol))
    for (ludyppvname in ludyppvnames) {
      subrast = subset(yrast,ludyppvname)
      plot(ecdf(values(subrast)), 
           xlab = paste0("Yield effect (%)"),
           ylab = "CDF", main = ludyppvname,
           cex = 0)
      
      zcolors = c("blue","green3","red")
      # txt = "Min. values\n"
      si = -5
      title(sub=paste0("Min. values (N. <=-",(cap*100.0),")"), adj=1, line=si, font=2)
      for (i in 1:length(zonecodes)) {
        zonecode = zonecodes[i]
        zcolor = zcolors[i]
        masked = mask(subrast,zones, maskvalue = zonecode, inverse = TRUE)
        lines(ecdf(values(masked)), 
              cex = 0, col = zcolor)
        minval = sprintf("%.1f",min(values(masked), na.rm = T))
        nval = sum(!is.na(values(masked)), na.rm = T)
        nbad = sum(values(masked) <= -(cap*100.0), na.rm = T)
        title(sub=paste0(names(zonecode)," : ",minval, "(",nbad,")"), adj=1, line=si+i, font=2, col.sub = zcolor)
      }
      # title(sub=txt, adj=1, line=-2, font=2)
    }
    tit = paste0(crop, " | ", scen, " ", year )
    title(tit, line = -1, outer = TRUE)
    par(mfrow = c(1,1))
    dev.off()
  }
  
  
  # FIXME This leads to a ~100dpi png, but elements dont scale well if we just set 'res' in png()
  png(filename = paste0(outfpref,crop,".png"), width = 1000, height = 500, unit = "px", pointsize = 16)
  print(plotobj)
  dev.off()
  png(filename = paste0(outfpref,crop,"_MARGINAL.png"), width = 1000, height = 500, unit = "px", pointsize = 16)
  print(plotobjmarg)
  dev.off()
  # png(filename = paste0(outfpref,crop,"_GDD.png"), width = 1500, height = 500*(length(dlus)-2), unit = "px", pointsize = 16)
  # print(tmap_arrange(plotobj,plotdgdd,plotdedd, ncol = 3))
  # dev.off()
  png(filename = paste0(outfpref,crop,"_VPD.png"), width = 1500, height = 500*(length(dlus)-2), unit = "px", pointsize = 16)
  print(tmap_arrange(plotobj,plotnoscale,plotanom, ncol = 3))
  dev.off()
  png(filename = paste0(outfpref,crop,"_ALL.png"), width = 3000, height = 500*(length(dlus)-2), unit = "px", pointsize = 16)
  print(tmap_arrange(plotobj,plotnoscale,plotdgdd,plotdtmp,plotdedd,plotdtmx, ncol = 6))
  dev.off()
  
}

# Write table
write.csv(outdata,paste0(wrtfolder,"/sr_",scen,"_",year,".csv"))

# breaks = classIntervals(na.omit(values(subrast)),n = 9, style="pretty")$brks
# plotdgdd <- tm_shape(subrast, bbox = bbox) +
#   tm_raster(palette = "YlOrRd", breaks = breaks,
#             title = expression(paste(Delta,"GDD"))) +
#   tm_shape(reg) + tm_borders() +
#   tm_legend(legend.text.size = 1.0, legend.title.size = 1.5,
#             legend.outside = T)
# plotdgdd
