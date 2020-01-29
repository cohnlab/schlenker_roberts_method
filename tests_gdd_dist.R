library(tidyverse)

infolder  = "xavier_dfs/"

crop = "Soybeans"

infname = paste0(infolder,crop,".climdata.rds")
data = readRDS(infname)
data$yvar = data$gdd30
# data$yvar = data$gdd30/data$ndays
# data$yvar = log(1+data$gdd30)/data$ndays
# data$yvar = log(1+data$gdd30/data$ndays)
data$xvar = data$tmaxmean

# fit = lm(yvar ~ xvar, data = data)
fit = nls(gdd30 ~ 0.001 + exp(b*tmaxmean), data = data)

summary(fit)
data$pred = predict(fit,data)

ggplot(data, aes(x = tmaxmean, y = gdd30, fill = ndays)) +
  geom_point(shape = 21, colour = "#00000000") +
  scale_fill_gradientn(colors = rainbow(5)) +
  theme_classic() +
  geom_line(aes(x = tmaxmean, y = pred, color = "black")) 


# ggplot(data, aes(x = xvar, y = yvar, fill = ndays)) +
#   geom_point(shape = 21, colour = "#00000000") +
#   scale_fill_gradientn(colors = rainbow(5)) +
#   theme_classic() +
#   geom_line(aes(x = xvar, y = pred, color = "black")) 

# plot = ggplot(data, aes(x = tmaxmean, y = gdd30/ndays, color = ndays)) +
#   geom_point() +
#   scale_color_gradientn(colors = rainbow(5)) +
#   theme_classic()
# plot
# plot + geom_line(aes(x = tmaxmean, y = pred$pred)) 
# plot + geom_line(data = pred, aes(x = tmaxmean, y = pred, color = "red")) 

# ggplot(data, aes(x = tmaxmean, y = gdd30, fill = ndays)) +
#   geom_point(shape = 21, colour = "#00000000") +
#   scale_fill_gradientn(colors = rainbow(5)) +
#   theme_classic() +
#   geom_line(aes(x = tmaxmean, y = pred$pred, color = "black")) 

