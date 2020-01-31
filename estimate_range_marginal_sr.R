# This reads a set of coefficent tables and generates estimates of yield impact for various LU change levels
Packages <- c("dplyr","tidyverse","data.table","sf","ncdf4","raster","fasterize","tmap","ggspatial","RColorBrewer")
lapply(Packages, library, character.only = TRUE)

tempbasefname = "../AgroServYield_coefficients/Inputs_coef/Tbaseline_tmn.csv"
tmaxbasefname = "../AgroServYield_coefficients/Inputs_coef/Tbaseline_tmx.csv"

betasfolder = "gdd_betas/"

# Points shapefile, reference raster and region of interest shapefile
shpfname = "GIS/COLROW30.shp"
reffname = "sheffield_calendars/Maize.crop.calendar.fill.nc"
regfname = "GIS/EstadosBR_IBGE_LLWGS84.shp"

# Path for output plots
outfolder = "plots/"

crops = c("Maize","Soybeans","Cotton")

# Number of days in the growing season used in SR20090
ndays = c(180,180,210)
names(ndays) <- crops

# Sensitivities of log yields to GDD and EDD. Found in digitizing/traced.ods
eddsens = c(-0.006435774107513,-0.005894686874659,-0.006753494536136)
gddsens = c(0.000317954802039,0.000403658172117,0.000907195445011)
names(eddsens) <- crops
names(gddsens) <- crops

# Effect of 1pp of LU change in maximum and mean temperatures
luefftemp = 0.91
luefftmax = 1.82
  
# FIXME do for a single year for now
year = 2030

# Get a scenario string from the filename
# scen = tools::file_path_sans_ext(basename(dttfname))

# Output plots
# outfpref = paste0(outfolder,"/estimate_range_",scen,"_",year)

# The values of change in land use to evaluate. Fractional
# dlus = seq(0,1,0.1)
dlus = c(0.00,0.05,0.10,0.20,0.3,0.5,1.0)

# Read baseline temp data
tempbase <- read.csv(tempbasefname) %>% dplyr::select(-X)
tmaxbase = read.csv(tmaxbasefname) %>% dplyr::select(-X)

# Make it long in respect to crops and join them
tempbase = gather(tempbase,"Crop","tempbase",-ID)
tmaxbase = gather(tmaxbase,"Crop","tmaxbase",-ID)
alldata <- left_join(tempbase,tmaxbase,by=c("ID","Crop"))

#FIXME: Just make a placeholder for GCC changes in temp and tmax
alldata$gccdtemp = 0.0
alldata$gccdtmax = 0.0


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

# # Computes the agroserv equation using columns in a data frame
# compute_agroserv <- function(df,dlu) {
#   with(alldata,
#        beta_1*(deltaT + agroservT*dlu) + 
#          beta_2*(deltaT + agroservT*dlu)^2 + 
#          beta_3*Tbaseline*(deltaT + agroservT*dlu) + 
#          beta_4*Tbaseline*(deltaT + agroservT*dlu)^2 + 
#          beta_5*AD*(deltaT + agroservT*dlu))
# }
# # Compute each level as a variable in the format dluX.XX
# dluvnames = paste0("dlu",sprintf("%.2f",dlus))
# for (i in 1:length(dlus)) {
#   alldata[dluvnames[i]] = compute_agroserv(df,dlus[i])
# }

ludtempvnames = paste0("ludtemp",sprintf("%.2f",dlus))
ludtmaxvnames = paste0("ludtmax",sprintf("%.2f",dlus))
ludgddvnames = paste0("ludgdd",sprintf("%.2f",dlus))
ludeddvnames = paste0("ludedd",sprintf("%.2f",dlus))
for (i in 1:length(dlus)) {
  alldata[ludtempvnames[i]] = luefftemp*dlus[i]
  alldata[ludtmaxvnames[i]] = luefftmax*dlus[i]
  
  alldata[ludgddvnames[i]] = eval_nxdd(betasgdd,(alldata["tempbase"]+alldata["gccdtemp"]+alldata[ludtempvnames[i]]) ) -
    eval_nxdd(betasgdd,alldata["tempbase"])
    # alldata$gccdtemp + alldata[ludtempvnames[i]]
}



# Read the shapefile and the reference raster
shp = st_read(shpfname)
ref = raster(reffname)

# FIXME Loop this
crop = "Maize"
# for (crop in levels(alldata$Crop)) {

# Read the coefficients for the T-GDD functions
betasedd = read.csv(paste0(betasfolder,crop,".betas.nEDD.csv"))
betasgdd = read.csv(paste0(betasfolder,crop,".betas.nGDD.csv"))

# Filter baselines for the crop of interest
cropbasedata <-basedata %>% filter(Crop ==  crop)

# FIXME Filter the year, we should be able to loop this too
# ydata <- cropdata %>% filter(ScenYear == year)

# Join with the shapefile
yshp = left_join(shp,ydata, by = c("COLROW30" = "ID"))

yrast <- rasterize(yshp, ref, field = dluvnames, fun = "last", background = NA_real_,
                   by = NULL)
# writeRaster(yrast,"poi.nc")

# Read the relevant overlay shapefile
reg = st_read(regfname)

# levs = c("dlu0.10","dlu0.20","dlu1.00")
# levs = c("dlu0.05","dlu0.10","dlu0.20","dlu1.00")
# subrast = subset(yrast,levs)
subrast = yrast
levs = names(subrast)

qmin = round(quantile(as.array(subrast),0.1, na.rm=T))
# breaks = seq(-5,5,0.5)
breaks = seq(qmin,-1*qmin,0.5)
pal = brewer.pal(n = length(breaks), name = "RdBu")

tit = paste0(crop, " | ", year, " in ", scen ," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))

plotobj <- tm_shape(subrast, bbox = reg) + 
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Yield (%)"))) + 
  tm_shape(reg) + tm_borders() +
  tm_legend(legend.text.size = 1.0,
            legend.outside = T) + 
  tm_layout(panel.show = T, panel.labels = levs,
            panel.label.size = 1.2,
            main.title.position = "center", main.title = tit)

# Repeating for the marginal effects
subrast = yrast - yrast$dlu0.00
names(subrast) <- names(yrast)
levs = names(subrast)

qmin = round(quantile(as.array(subrast),0.1, na.rm=T))
# breaks = seq(-5,5,0.5)
breaks = seq(qmin,-1*qmin,0.5)
pal = brewer.pal(n = length(breaks), name = "RdBu")

tit = paste0(crop, " (Marginal: X-dlu0.00) | ", year, " in ", scen ," | Min: ", sprintf("%.2f",min(as.array(subrast), na.rm=T)))

plotobjmarg <- tm_shape(subrast, bbox = reg) + 
  tm_raster(palette = pal, breaks = breaks,
            title = expression(paste(Delta,"Yield (%)"))) + 
  tm_shape(reg) + tm_borders() +
  tm_legend(legend.text.size = 1.0,
            legend.outside = T) + 
  tm_layout(panel.show = T, panel.labels = levs,
            panel.label.size = 1.2,
            main.title.position = "center", main.title = tit)


# FIXME This leads to a ~100dpi png, but elements dont scale well if we just set 'res' in png()
png(filename = paste0(outfpref,crop,".png"), width = 800, height = 500, unit = "px", pointsize = 16)
print(plotobj)
dev.off()
png(filename = paste0(outfpref,crop,"_MARGINAL.png"), width = 800, height = 500, unit = "px", pointsize = 16)
print(plotobjmarg)
dev.off()
# }

# spplot(subset(yrast,c("dlu0.10","dlu0.20")), 
#      xlim = c(-80,-30), ylim = c(-33.0, 10.0),
#      )

# tm_shape(subrast, bbox = reg) + 
#   tm_raster(palette = "RdBu", n = 9, midpoint = 0.0) + 
#   tm_shape(reg) + tm_borders()
# plotobj <- tm_shape(subrast) + 
#   tm_raster(palette = pal, breaks = breaks,
#             title = expression(paste(Delta,"Yield (%)"))) +
#   tm_legend(legend.text.size = 1.0,
#             legend.outside = T) + 
#   tm_layout(panel.show = T, panel.labels = levs,
#             panel.label.size = 1.2,
#             main.title.position = "center", main.title = tit)
