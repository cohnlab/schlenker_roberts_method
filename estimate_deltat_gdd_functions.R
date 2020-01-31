library(tidyverse)
library(ggplot2)
library(gridExtra)
library(lspline)

infolder  = "sheffield_dfs/"
outfolder = "gdd_betas/"

# years = 1991:1993
years = 1991:2008

crops = c("Maize","Soybeans","Cotton")

# Each crop will use different tresholds for GDD and EDD.
# Give lower limit EDD names for GDD, and later GDD can be calculated as EDDl-EDDu
# Values taken from Schlenker and Roberts can be found in digitized/traced.ods
eddvnames = c("edd29","edd30","edd32")
gddvnames = c("edd10","edd10","edd15")
names(eddvnames) <- crops
names(gddvnames) <- crops

# Cuts (hinges, breaks) of the spline approximation
cuts = seq(9,37,2)

# Functions

#FIXME: This function do not allow for a decreasing last piece to avoid problems with rank-defficiency
lspline_to_betas <- function(fit,cuts) {
  vname = all.vars(formula(fit))[2]
  betas = data.frame(
    sta = -999,
    en = cuts[1],
    intr = fit$coefficients["(Intercept)"],
    slp = fit$coefficients[2]
  )
  cuts = c(cuts,999) # Add a number to be the last cut
  for (i in 2:(length(cuts))) {
    sta = cuts[i-1]
    en = cuts[i]
    slp = fit$coefficients[i+1]
    p = data.frame(sta)
    names(p) <- vname
    intr = predict(fit,p) - slp*sta
    if (is.na(intr) | (slp<0.0)) { # FIXME: We force the last one to be positive
      slp = betas$slp[i-1]
      intr = betas$intr[i-1]
    }
    betas = rbind(betas,
                  data.frame(sta=sta,en=en,intr=intr,slp=slp))
    
  }
  # betas = rbind(betas,
  # data.frame(sta=en,en=999,intr=intr,slp=slp))
  rownames(betas) <- NULL
  return(betas)
}
add_lines_betas <- function(betas) {
  for (i in 1:nrow(betas)) {
    pdata = seq(betas$sta[i],betas$en[i],0.1)
    lines(pdata,betas$intr[i] + betas$slp[i]*pdata,col='red',lwd=2.0)
  }
}

par(mfrow = c(length(crops),2))
# crop = crops[1]
for (crop in crops) {
  infname = paste0(infolder,crop,".climdata.rds")
  data = readRDS(infname) %>% filter(year %in% years)
  
  data["EDD"] <- data[eddvnames[crop]]
  data["GDD"] <- data[gddvnames[crop]] - data["EDD"]
  data["nEDD"] <- data["EDD"]/data["ndays"]
  data["nGDD"] <- data["GDD"]/data["ndays"]
  
  
  # Fit nEDD to Tmax and nGDD to Tavg
  fitedd = lm(nEDD ~ lspline(tmaxmean,cuts), data = data)
  fitgdd = lm(nGDD ~ lspline(tempmean,cuts), data = data)
  # Convert to table of intercepts and slopes using the function lspline_to_betas
  betasedd = lspline_to_betas(fitedd,cuts)
  betasgdd = lspline_to_betas(fitgdd,cuts)
  
  # Add the plots and lines
  with(data,plot(tmaxmean,nEDD,main=crop,pch='.', cex = 2.0))
  add_lines_betas(betasedd)
  with(data,plot(tempmean,nGDD,main=crop,pch='.', cex = 2.0))
  add_lines_betas(betasgdd)
  
  # Create output folder
  dir.create(file.path(outfolder), showWarnings = FALSE)
  
  # Write each crop and nEDD/nGDD piecewise coefficients and hinges
  write.csv(betasedd, paste0(outfolder, crop, ".betas.nEDD.csv"), row.names = F)
  write.csv(betasgdd, paste0(outfolder, crop, ".betas.nGDD.csv"), row.names = F)
}

#' b0 = "$\\beta_0 + \\frac{1}{2}$"
#' 
#' #'
#' {{b0}}