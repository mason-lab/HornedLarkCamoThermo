#The following script was used to analyze the relationship between cooling benefits and the temperature, aridity, and seasonality PCs

#import libraries
library(minpack.lm)
library(lme4)
library(MuMIn)

#import data
data <- read.csv("~/Data/Endoscape_cooling_benefit_data.csv")


#model selection
fit <- nls(Delta_cooling_avg_dorsal~a*exp(b*SeasonalityPC),data=data,start=list(a=0.07,b=0.000054))
fit2 <- nls(Delta_cooling_avg_dorsal~a*exp(b*AridityPC),data=data,start=list(a=0.04,b=0.0027))
fit3 <- nls(Delta_cooling_avg_dorsal~a*exp(b*TemperaturePC),data=data,start=list(a=0.05,b=0.006))
fit4 <- lm(Delta_cooling_avg_dorsal~1,data=data)
model.sel(fit,fit2,fit3,fit4)

summary(fit2)

#Bootstrapping method to generate 95% confidence intervals

#bootstrapping method: seasonal
x <- data$SeasonalityPC
y <- data$Delta_cooling_avg_dorsal
new.x <- seq.int(-4497, 7664, by = 176)

estimate <- function(ind){
  x <- x[ind]
  y <- y[ind]
  m1 <- nls(y ~ a*exp(b*x), start=list(a=0.07,b=0.000054),
            control = nls.control(maxiter = 500, warnOnly = TRUE))
  predict(m1, newdata = list(x = new.x))
}
predict0 <- predict(fit, newdata = list(time = all.time))
predict1 <- replicate(1000, estimate(sample.int(70, replace = TRUE)))
intervals <- apply(predict1, 1, quantile, probs = c(0.025, 0.975))
final <- rbind(predict0, intervals,new.x)
write.csv(x = final,file='Output/seasonal.csv')

#bootstrapping method: aridity
x <- data$AridityPC
y <- data$Delta_cooling_avg_dorsal
new.x <- seq.int(-1128, 480, by = 23)

estimate <- function(ind){
  x <- x[ind]
  y <- y[ind]
  m1 <- nls(y ~ a*exp(b*x), start=list(a=0.04,b=0.0027),
            control = nls.control(maxiter = 500, warnOnly = TRUE))
  predict(m1, newdata = list(x = new.x))
}
predict0 <- predict(fit2, newdata = list(time = all.time))
predict1 <- replicate(1000, estimate(sample.int(70, replace = TRUE)))
intervals <- apply(predict1, 1, quantile, probs = c(0.025, 0.975))
final <- rbind(predict0, intervals,new.x)
write.csv(x = final,file='Output/aridity.csv')

#bootstrapping method: temperature
x <- data$TemperaturePC
y <- data$Delta_cooling_avg_dorsal
new.x <- seq.int(-383, 258, by = 9.2)

estimate <- function(ind){
  x <- x[ind]
  y <- y[ind]
  m1 <- nls(y ~ a*exp(b*x), start=list(a=0.05,b=0.006),
            control = nls.control(maxiter = 500, warnOnly = TRUE))
  predict(m1, newdata = list(x = new.x))
}
predict0 <- predict(fit3, newdata = list(time = all.time))
predict1 <- replicate(1000, estimate(sample.int(70, replace = TRUE)))
intervals <- apply(predict1, 1, quantile, probs = c(0.025, 0.975))
final <- rbind(predict0, intervals,new.x)
write.csv(x = final,file='Output/temperature.csv')
