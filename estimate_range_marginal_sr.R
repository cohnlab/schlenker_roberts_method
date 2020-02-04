# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer")
lapply(Packages, library, character.only = TRUE)

tempbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmn.csv"
tmaxbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmx.csv"

betasfolder = "gdd_betas/"

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
  msklim = 0.001
} else {
  regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"
  msklim = 0.000
}


brmskfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"


# Path for output plots
if (iglobal) {
  outfolder = "plots_sr_memo_global/"
} else {
  outfolder = "plots_sr_memo/"
}

crops = c("Maize","Soybeans","Cotton")
cropstrings = c("maize","soybean","cotton")
names(cropstrings) <- crops

# Number of days in the growing season used in SR20090
# ndays = c(180,180,210)
ndays = c(120,120,120)
names(ndays) <- crops

# Cap for logY impacts, both positive and negative
cap = 0.5

# Sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Effect of 1pp of LU change in maximum and mean temperatures
luefftemp = 1.57*(0.5/0.7)
luefftmax = 3.15*(0.5/0.7)

# FIXME: Fixed backgound warming.Substitute for a RCP file
bgtemp = 1.0
bgtmax = 1.0
basetit = paste0("Background ",round(bgtemp),"\u00B0C warming")


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
dlus = c(0.0,0.3)
# dlus = c(.00,0.05,0.10,0.20,0.3,0.5,1.0)
# dlus = c(0.00,0.05,0.10,0.20,0.3,1.0)

# Grid for arranging the plots
gnrow = 1
gncol = length(dlus)-1

# Maximum number of bins in plot. More than that, it increases bin size
maxbreaks = 30

# Override quantile based breaks anyway
ioverbreaks = TRUE
overbreaks = seq(-30,30,5)

# Create output folder
dir.create(outfolder, showWarnings = F)

# Read baseline temp data
tempbase <- read.csv(tempbasefname) %>% dplyr::select(-X)
tmaxbase = read.csv(tmaxbasefname) %>% dplyr::select(-X)

# Make it long in respect to crops and join them
tempbase = gather(tempbase,"Crop","tempbase",-ID)
tmaxbase = gather(tmaxbase,"Crop","tmaxbase",-ID)
alldata <- left_join(tempbase,tmaxbase,by=c("ID","Crop"))

#FIXME: Just make a placeholder for GCC changes in temp and tmax
alldata$gccdtemp = bgtemp
alldata$gccdtmax = bgtmax


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
  bbox = globalbbox
} else {
  bbox = reg
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
  
  # Here we actually apply the equations, setting variable names prepended with the different LU levels
  ludtempvnames = paste0("ludtemp",sprintf("%.2f",dlus))
  ludtmaxvnames = paste0("ludtmax",sprintf("%.2f",dlus))
  ludgddvnames = paste0("dgddDLU",sprintf("%.2f",dlus))
  ludeddvnames = paste0("deddDLU",sprintf("%.2f",dlus))
  ludlogyvnames = paste0("dlogyDLU",sprintf("%.2f",dlus))
  ludyppvnames = paste0("dyppDLU",sprintf("%.2f",dlus))
  for (i in 1:length(dlus)) {
    # Changes in temp, tmax
    ydata[ludtempvnames[i]] = luefftemp*dlus[i]
    ydata[ludtmaxvnames[i]] = luefftmax*dlus[i]
    
    # Changes in GDD, EDD, multiplying by ndays
    ydata[ludgddvnames[i]] = (eval_nxdd(betasgdd,(ydata["tempbase"]+ydata["gccdtemp"]+ydata[ludtempvnames[i]]) ) -
                                eval_nxdd(betasgdd,ydata["tempbase"]))*ndays[crop]
    ydata[ludeddvnames[i]] = (eval_nxdd(betasedd,(ydata["tmaxbase"]+ydata["gccdtmax"]+ydata[ludtmaxvnames[i]]) ) -
                                eval_nxdd(betasedd,ydata["tmaxbase"]))*ndays[crop]
    
    # Apply the GDD and EDD sensitivities
    ydata[ludlogyvnames[i]] = gddsens[crop]*ydata[ludgddvnames[i]] +
      eddsens[crop]*ydata[ludeddvnames[i]]
    
    #FIXME: Cap impacts at cap
    ydata[[ludlogyvnames[i]]][ydata[[ludlogyvnames[i]]] > cap] = cap
    ydata[[ludlogyvnames[i]]][ydata[[ludlogyvnames[i]]] < -cap] = -cap
    
    # Get a percent value
    ydata[ludyppvnames[i]] = ydata[ludlogyvnames[i]]*100.0
  }
  
  # Join with the shapefile
  yshp = left_join(shp,ydata, by = c("COLROW30" = "ID"))
  
  yrast <- rasterize(yshp, ref, field = ludyppvnames, fun = "last", background = NA_real_,
                     by = NULL)
  if (imask) {
    yrast = mask(yrast,msk)
  }
  
  # levs = c("dlu0.10","dlu0.20","dlu1.00")
  # levs = c("dlu0.05","dlu0.10","dlu0.20","dlu1.00")
  # subrast = subset(yrast,levs)
  subrast = subset(yrast,names(yrast)[-1])
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
  
  # tit = paste0(crop, " | ", year, " in ", scen ," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  tit = paste0(crop, " | ", basetit," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  
  plotobj <- tm_shape(subrast, bbox = bbox) + 
    tm_raster(palette = pal, breaks = breaks,
              title = expression(paste(Delta,"Yield (%)"))) + 
    tm_shape(reg) + tm_borders() +
    tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
              legend.outside = T) + 
    tm_layout(panel.show = F, panel.labels = paste0("dlu",dlus),
              panel.label.size = 1.2,
              main.title.position = "center", main.title = tit, main.title.size = 1.0) +
    tm_facets(nrow = gnrow, ncol = gncol)
  
  # Repeating for the marginal effects
  subrast = yrast - yrast$dyppDLU0.00
  names(subrast) <- names(yrast)
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
  
  # tit = paste0(crop, " (Marginal: X-dlu0.00) | ", year, " in ", scen ," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  #tit = paste0(crop, " (Marginal) | ", basetit ," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  if (iglobal) {
    tit = paste0("Impact to ", cropstrings[crop], " yield of 30 percentage point loss of native vegetation")
    
  } else { 
    tit = paste0("Impact to ", cropstrings[crop], " yield of 30 percentage point loss of native vegetation - Brazil")
    
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
  
  
  
  # FIXME This leads to a ~100dpi png, but elements dont scale well if we just set 'res' in png()
  png(filename = paste0(outfpref,crop,".png"), width = 1000, height = 500, unit = "px", pointsize = 16)
  print(plotobj)
  dev.off()
  png(filename = paste0(outfpref,crop,"_MARGINAL.png"), width = 1000, height = 500, unit = "px", pointsize = 16)
  print(plotobjmarg)
  dev.off()
  
}
