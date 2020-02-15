# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("Rcpp","dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer")
lapply(Packages, library, character.only = TRUE)

# Land use change values for the polygon grid
dlufname = "C:/Users/eriad/Downloads/newAGT1_DeltaLUCs.csv"

# File with equivalencies for conversion from the polygon grid to the points grid
convfname = "../AgroServYield_coefficients/Inputs_coef/grids.csv"

tempbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmn.csv"
tmaxbasefname = "../AgroServYield_coefficients/Inputs_coef/CRU_Sacks_Zarc/Tbaseline_tmx.csv"

betasfolder = "gdd_betas/"

# Points shapefile, reference raster and region of interest shapefile
shpfname = "GIS/COLROW30.shp"
reffname = "../AgroServYield_coefficients/Inputs_Tbaseline/cru_tmp.nc"

iglobal = FALSE
if (iglobal) {
  regfname = "GIS/bbox_world.shp"
} else {
  regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"
}


# Path for output plots
outfolder = "plots_sr_memo/"

crops = c("Maize","Soybeans","Cotton")
# crops = c("Soybeans")

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
bgtemp = 0.0
bgtmax = 0.0
basetit = paste0("Background ",round(bgtemp),"\u00B0C warming")

# Get a scenario string from the filename
# scen = tools::file_path_sans_ext(basename(dttfname))
scen = paste0("deltaT",bgtemp)

# Output plots
outfpref = paste0(outfolder,"/globiom_lu_",scen)

# The values of change in land use to evaluate. Fractional
# dlus = c(0.00,0.05,0.10,0.20,0.3,1.0)
years = c(2050)

# Grid for arranging the plots
gnrow = 1
gncol = length(years)

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

# Reading LUC data and using the conversion table to give it polygon IDs
convdata = read.csv(convfname)
convdata$COLROW30 <- as.character(convdata$final.COLROW30) # Improves compatibility

dludata = read.csv(dlufname)
dludata <- left_join(dludata,convdata,by = c( "ID" = "final.ID"))

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
  outvals = cpp_eval_nxdd(as.matrix(betas),as.vector(tvec))
  # outvals[is.na(tvec)] <- NA
  return(outvals)
}

# Read the reference raster
ref = raster(reffname)

# Read the shapefile, join it with LUC values and rasterize them
shp = st_read(shpfname)
shp$COLROW30 <- as.character(shp$COLROW30) # Improves compatibility

shp = left_join(shp,dludata)

lurast <- rasterize(shp, ref, field = grep("X[0-9]{4}",names(shp), value = T), fun = mean, background = NA_real_,
                  by = NULL)

# Read the relevant overlay shapefile
reg = st_read(regfname)

# LUC plot plot
luplot = tm_shape(subset(lurast,paste0("X",years)), bbox = reg) + 
  tm_raster(breaks = seq(0,100,10),
            interval.closure = "left",
            title = expression(paste(Delta,"LU (pp)"))) + 
  tm_shape(reg) + tm_borders() +
  tm_legend(legend.text.size = 1.0, legend.title.size = 1.5, 
            legend.outside = T) + 
  tm_layout(panel.show = T, panel.labels = years,
            legend.format = list("text.less.than" = "<"),
            panel.label.size = 1.2,
            main.title.position = "center", main.title = "Percentage point loss of\n native vegetation relative to size of pixel",
            main.title.size = 1.0) #+
  #tm_facets(nrow = gnrow, ncol = gncol)
luplot




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
  
  # Add dlu values to ydata
  ydata = left_join(ydata,dludata,by = c("ID" = "COLROW30"))
  
  # Here we actually apply the equations, setting variable names prepended with the different LU levels
  dlus = years
  luvnames = paste0("X",dlus)
  names(luvnames) = dlus
  ludtempvnames = paste0("ludtemp",sprintf("%4i",dlus))
  ludtmaxvnames = paste0("ludtmax",sprintf("%4i",dlus))
  ludgddvnames = paste0("dgddDLU",sprintf("%4i",dlus))
  ludeddvnames = paste0("deddDLU",sprintf("%4i",dlus))
  ludlogyvnames = paste0("dlogyDLU",sprintf("%4i",dlus))
  ludyppvnames = paste0("dyppDLU",sprintf("%4i",dlus))
  for (i in 1:length(dlus)) {
    # Changes in temp, tmax
    ydata[ludtempvnames[i]] = luefftemp*ydata[luvnames[i]]*0.01
    ydata[ludtmaxvnames[i]] = luefftmax*ydata[luvnames[i]]*0.01
    
    # Changes in GDD, EDD, multiplying by ndays
    ydata[ludgddvnames[i]] = (eval_nxdd(betasgdd,(ydata["tempbase"]+ydata["gccdtemp"]+ydata[ludtempvnames[i]]) ) -
                                eval_nxdd(betasgdd,ydata["tempbase"]))*ndays[crop]
    ydata[ludeddvnames[i]] = (eval_nxdd(betasedd,(ydata["tmaxbase"]+ydata["gccdtmax"]+ydata[ludtmaxvnames[i]]) ) -
                                eval_nxdd(betasedd,ydata["tmaxbase"]))*ndays[crop]
    
    #FIXME: Filter out NAs here since eval_nxdd interprets them as zeros. 
    #Should be done in the wrapper but something breaks if implemented simply.
    ydata[[ludgddvnames[i]]][is.na(ydata[ludtempvnames[i]])] <- NA
    ydata[[ludeddvnames[i]]][is.na(ydata[ludtempvnames[i]])] <- NA
    
    
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
  
  yrast <- rasterize(yshp, ref, field = ludyppvnames, fun = mean, background = NA_real_,
                     by = NULL)
  
  subrast = yrast
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
  tit = paste0(crop," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))
  tit = paste0("Impact of natural vegetation loss on\n",crop," yields")
  
  plotobj <- tm_shape(subrast, bbox = reg) + 
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
  png(filename = paste0(outfpref,crop,".png"), width = 600, height = 300, unit = "px", pointsize = 16)
  print(plotobj)
  dev.off()
  
}

png(filename = paste0(outfpref,"LUC",".png"), width = 600, height = 300, unit = "px", pointsize = 16)
print(luplot)
dev.off()