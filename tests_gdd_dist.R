library(tidyverse)
library(earth)
library(segmented)
library(gridExtra)
library(lspline)
library(splines)

# infolder  = "xavier_dfs/"
infolder  = "sheffield_dfs/"

crop = "Maize"

# Option to filter years for tinkering
# years = 2000:2002
years = 1991:1993

infname = paste0(infolder,crop,".climdata.rds")
data = readRDS(infname) %>% filter(year %in% years)
data$nedd30 = data$edd30/data$ndays
data$ngdd1030 = data$gdd1030/data$ndays

# # data$yvar = data$edd30/data$ndays
# # data$yvar = log(1+data$edd30)/data$ndays
# # data$yvar = log(1+data$edd30/data$ndays)
# data$xvar = data$tempmean
# 
# # fit = lm(yvar ~ xvar, data = data)
# # fit = nls(nedd30 ~ 0.000001 + exp(b*tempmean),
# #           start = c(b=0.000001),
# #           data = data)
# 
# # fit = earth(nedd30 ~ tempmean,
# #             degree = 2,nk =3,
# #             # allowed = function(degree, pred, parents, namesx) namesx[pred] != ,
# #             data = data)
# # fit = segmented(lm(nedd30 ~ tempmean + I(tempmean^2), data = data),
# #                 segZ = ~ tempmean,
# #                 npsi = 1,
# #                 # psi = list(tempmean = c(25,30)),
# #             data = data)
# fit = earth(nedd30 ~ tempmean,
#             # degree = 2,nk =3,
#             # Auto.linpreds = FALSE,
#             # thresh=0, penalty=-1,
#             # allowed = function(degree, pred, parents, namesx) namesx[pred] != ,
#             data = data)


fit = segmented(lm(nedd30 ~tempmean, data = data),
                segZ = ~tempmean,
                psi = list(tempmean = c(20)))

cuts = seq(9,37,2)

fit = lm(nedd30 ~ lspline(tmaxmean,cuts), data = data)
with(data,plot(tmaxmean,nedd30,pch='.', cex = 2.0))
pdata = data.frame(tmaxmean = seq(0,50,0.1))
pdata$pred = predict(fit,pdata)
with(pdata,lines(tmaxmean,pred, col="orange",lwd=1.2))

summary(fit)
fit$coefficients

cuts = seq(9,37,2)

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
    if (is.na(intr)) {
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
lspline_to_betas(fit,cuts)

betas = lspline_to_betas(fit,cuts)

add_full_line <- function(betas,i) {
  pdata = seq(0,50,0.1)
  lines(pdata,betas$intr[i] + betas$slp[i]*pdata,col='red',lwd=2.5)
}

add_part_line <- function(betas,i) {
  pdata = seq(betas$sta[i],betas$en[i],0.1)
  lines(pdata,betas$intr[i] + betas$slp[i]*pdata,col='green',lwd=2.0)
}


with(data,plot(tmaxmean,nedd30,pch='.', cex = 2.0))
for (i in 1:nrow(betas)) {
  add_full_line(betas,i)
}
for (i in 1:nrow(betas)) {
  add_part_line(betas,i)
}
  
cuts = seq(9,37,2)

fit = lm(ngdd1030 ~ lspline(tempmean,cuts), data = data)
with(data,plot(tempmean,ngdd1030,pch='.', cex = 2.0))
pdata = data.frame(tempmean = seq(0,50,0.1))
pdata$pred = predict(fit,pdata)
with(pdata,lines(tempmean,pred, col="red"))

betas = lspline_to_betas(fit,cuts)

summary(fit)
fit$coefficients



part = filter(data, tempmean < 20)
fit = lm(nedd30 ~ poly(tempmean,1), data = part)
part$pred = predict(fit,part)
with(part,points(tempmean,pred, col="red"))


part = filter(data, tempmean >= 20)
fit = lm(nedd30 ~ poly(tempmean,2), data = part)
part$pred = predict(fit,part)
with(part,points(tempmean,pred, col="red"))

cuts = seq(5,35,2)
with(data,plot(tempmean,ngdd1030,pch='.', cex = 2.0))
for (i in  seq(1,length(cuts)-1)) {
  part = filter(data, tempmean >= cuts[i] & tempmean <= cuts[i+1])
  fit = lm(ngdd1030 ~ tempmean, data = part)
  pdata = data.frame(tempmean = seq(cuts[i],cuts[i+1],0.1))
  pdata$pred = predict(fit,pdata)
  with(pdata,lines(tempmean,pred, col="red"))
}

cuts = seq(5,40,2)
# with(data,plot(tmaxmean,nedd30,pch='.', cex = 2.0,xlim=c(20,35),ylim=c(0,1.5)))
with(data,plot(tmaxmean,nedd30,pch='.', cex = 2.0))
for (i in  seq(1,length(cuts)-1)) {
  part = filter(data, tmaxmean >= cuts[i] & tmaxmean <= cuts[i+1])
  fit = lm(nedd30 ~ tmaxmean, data = part)
  pdata = data.frame(tmaxmean = seq(cuts[i],cuts[i+1],0.1))
  pdata$pred = predict(fit,pdata)
  with(pdata,lines(tmaxmean,pred, col="red",lwd=2.0))
}


# 
# ggplot(data, aes(x = xvar, y = nedd30, fill = ndays)) +
#   geom_point(shape = 21, colour = "#00000000") +
#   scale_fill_gradientn(colors = rainbow(5)) +
#   theme_classic() +
#   geom_line(aes(x = xvar, y = pred, color = "black", size=1.2))


plot = ggplot(data, aes(x = tmaxmean, y = nedd30, fill = trngmean)) +
  geom_point(shape = 21, colour = "#00000000") +
  scale_fill_gradientn(colors = rainbow(5)) +
  theme_classic()
plot + facet_grid(vars(cut(trngmean,quantile(trngmean))))


ggplot(data, aes(x = tempmean, y = nedd30)) +
  geom_point() +
  theme_classic()


plot1 = ggplot(data, aes(x = tempmean, y = nedd30)) +
  geom_point() +
  theme_classic()
plot2 = ggplot(data, aes(x = tminmean, y = nedd30)) +
  geom_point() +
  theme_classic()
plot3 = ggplot(data, aes(x = tmaxmean, y = nedd30)) +
  geom_point() +
  theme_classic()
grid.arrange(plot1,plot2,plot3)


plot1 = ggplot(data, aes(x = tempmean, y = ngdd1030)) +
  geom_point() +
  theme_classic()
plot2 = ggplot(data, aes(x = tminmean, y = ngdd1030)) +
  geom_point() +
  theme_classic()
plot3 = ggplot(data, aes(x = tmaxmean, y = ngdd1030)) +
  geom_point() +
  theme_classic()
grid.arrange(plot1,plot2,plot3)

with(data, plot(tempmean,ngdd1030))
with(data, plot(tmaxmean,nedd30))

ggplot(data, aes(x = tmaxmean, y = edd30, fill = ndays)) +
  geom_point(shape = 21, colour = "#00000000") +
  scale_fill_gradientn(colors = rainbow(5)) +
  theme_classic()


# ggplot(data, aes(x = xvar, y = yvar, fill = ndays)) +
#   geom_point(shape = 21, colour = "#00000000") +
#   scale_fill_gradientn(colors = rainbow(5)) +
#   theme_classic() +
#   geom_line(aes(x = xvar, y = pred, color = "black")) 

# plot = ggplot(data, aes(x = xvar, y = edd30/ndays, color = ndays)) +
#   geom_point() +
#   scale_color_gradientn(colors = rainbow(5)) +
#   theme_classic()
# plot
# plot + geom_line(aes(x = xvar, y = pred$pred)) 
# plot + geom_line(data = pred, aes(x = xvar, y = pred, color = "red")) 

# ggplot(data, aes(x = xvar, y = edd30, fill = ndays)) +
#   geom_point(shape = 21, colour = "#00000000") +
#   scale_fill_gradientn(colors = rainbow(5)) +
#   theme_classic() +
#   geom_line(aes(x = xvar, y = pred$pred, color = "black")) 

