library(tidyverse)
library(vars)
library(readxl)
library(lubridate)
library(hrbrthemes)
library(stringr)
library(forecast)

month <- read_excel("my_uzb_data.xlsx", sheet = "monthly")
quarter <- read_excel("my_uzb_data.xlsx", sheet = "quarterly")

money <- month[,c("m0","m1","m2","m2x")]
money <- ts(money[seq(2,dim(money)[1],3),], start=c(2013,1), frequency=4)

my_plot <- function(data) {
  for (col in colnames(data)) {
    par(mfrow=c(2,1))
    plot(data[,col], main=col)
    acf(data[,col], na.action=na.pass, main="")
  }
}

m_ts <- ts(month[,c(2:5,7:8)], start = c(2013,2), frequency = 12)
q_ts <- ts(quarter[,c(2,4,5)], start = c(2010,1), frequency = 4)

q <- ts(aggregate.ts(m_ts, nfrequency = 4, FUN=mean), start = c(2013,1), frequency = 4)

end.vars <- ts.intersect(q_ts[,c("gdp_real","cpi")],q[,"lending_rate"],money)
end.vars <- na.omit(end.vars)
colnames(end.vars) <- c("gdp_real","cpi","lending_rate",colnames(money))
end.vars[,"cpi"] <- end.vars[,"cpi"] * 100

ex.vars <- ts.intersect(q_ts[,"remittances"],q[,-4])
ex.vars <- na.omit(ex.vars)
colnames(ex.vars) <- c("remittances",colnames(q[,-4]))

mape <- function(test, pred) return(mean(abs((test-pred)/test))*100)

# M0

train <- window(end.vars, end=c(2020,4))[,c("gdp_real","cpi","lending_rate","m0")]
test <- window(end.vars, start=c(2021,1))[,"m0"]
train_x <- window(ex.vars, start=c(2016,1), end=c(2020,4))[,"cotton"]
train_x <- matrix(train_x,length(train_x),1)
test_x <- window(ex.vars, start=c(2021,1), end=c(2021,3))[,"cotton"]
test_x <- matrix(test_x,length(test_x),1)

mod0 <- VAR(train, p=1, type="both", season=4, exogen=train_x)
summary(mod0)
pred0 <- predict(mod0, n.ahead=3, dumvar=test_x)$fcst[["m0"]][,c("fcst","lower","upper")]

data <- data.frame(
  quarter = as.Date("2016-02-15") + (0:22)*90,
  value = as.vector(end.vars[,"m0"])
)
data[c("forecast","lower","upper")] <- NA
data[(dim(data)[1]-2):dim(data)[1],c("forecast","lower","upper")] <- pred0
mape0 <- round(mape(test,pred0[,1]),3)

ggplot(data, aes(x=quarter)) +
  geom_line(aes(y=value)) +
  geom_line(aes(y=upper), color="red", linetype="dashed", size=1) +
  geom_line(aes(y=forecast), color="blue", linetype="longdash", size=1) +
  geom_line(aes(y=lower), color="red", linetype="dashed", size=1) +
  theme_light() + labs(x="", y="") + 
  ggtitle(str_glue("Наличная денежная масса M0 (MAPE = {mape0})"))

# plot(predict(mod0, n.ahead=3, dumvar=test_x), names="m0",
#      main=str_glue("Денежная база M0 (MAPE = {mape0})"))

# M1

train <- window(end.vars, end=c(2020,4))[,c("gdp_real","cpi","lending_rate","m1")]
test <- window(end.vars, start=c(2021,1))[,"m1"]
train_x <- window(ex.vars, start=c(2016,1), end=c(2020,4))[,"cotton"]
train_x <- matrix(train_x,length(train_x),1)
test_x <- window(ex.vars, start=c(2021,1), end=c(2021,3))[,"cotton"]
test_x <- matrix(test_x,length(test_x),1)

mod1 <- VAR(train, p=1, type="both", season=4, exogen=train_x)
summary(mod1, equations="m1")
pred1 <- predict(mod1, n.ahead=3, dumvar=test_x)$fcst[["m1"]][,c("fcst","lower","upper")]

data <- data.frame(
  quarter = as.Date("2016-02-15") + (0:22)*90,
  value = as.vector(end.vars[,"m1"])
)
data[c("forecast","lower","upper")] <- NA
data[(dim(data)[1]-2):dim(data)[1],c("forecast","lower","upper")] <- pred1
mape1 <- round(mape(test,pred1[,1]),3)

ggplot(data, aes(x=quarter)) +
  geom_line(aes(y=value)) +
  geom_line(aes(y=upper), color="red", linetype="dashed", size=1) +
  geom_line(aes(y=forecast), color="blue", linetype="longdash", size=1) +
  geom_line(aes(y=lower), color="red", linetype="dashed", size=1) +
  theme_light() + labs(x="", y="") + 
  ggtitle(str_glue("Узкая денежная масса M1 (MAPE = {mape1})"))

# M2

train <- window(end.vars, end=c(2020,4))[,c("gdp_real","cpi","lending_rate","m2")]
test <- window(end.vars, start=c(2021,1))[,"m2"]
train_x <- window(ex.vars, start=c(2016,1), end=c(2020,4))[,"cotton"]
train_x <- matrix(train_x,length(train_x),1)
test_x <- window(ex.vars, start=c(2021,1), end=c(2021,3))[,"cotton"]
test_x <- matrix(test_x,length(test_x),1)

mod2 <- VAR(train, p=1, type="both", season=4, exogen=train_x)
summary(mod2, equation="m2")
pred2 <- predict(mod2, n.ahead=3, dumvar=test_x)$fcst[["m2"]][,c("fcst","lower","upper")]

data <- data.frame(
  quarter = as.Date("2016-02-15") + (0:22)*90,
  value = as.vector(end.vars[,"m2"])
)
data[c("forecast","lower","upper")] <- NA
data[(dim(data)[1]-2):dim(data)[1],c("forecast","lower","upper")] <- pred2
mape2 <- round(mape(test,pred2[,1]),3)

ggplot(data, aes(x=quarter)) +
  geom_line(aes(y=value)) +
  geom_line(aes(y=upper), color="red", linetype="dashed", size=1) +
  geom_line(aes(y=forecast), color="blue", linetype="longdash", size=1) +
  geom_line(aes(y=lower), color="red", linetype="dashed", size=1) +
  theme_light() + labs(x="", y="") + 
  ggtitle(str_glue("Широкая денежная масса в национальном определении M2\n(MAPE = {mape2})"))

# M2X

train <- window(end.vars, end=c(2020,3))[,c("gdp_real","cpi","lending_rate","m2x")]
test <- window(end.vars, start=c(2020,4))[,"m2x"]
train_x <- window(ex.vars, start=c(2016,1), end=c(2020,3))[,"cotton"]
train_x <- matrix(train_x,length(train_x),1)
test_x <- window(ex.vars, start=c(2020,4), end=c(2021,3))[,"cotton"]
test_x <- matrix(test_x,length(test_x),1)

mod3 <- VAR(train, p=1, type="both", season=4, exogen=train_x)
summary(mod3, equation="m2x")
pred3 <- predict(mod3, n.ahead=4, dumvar=test_x)$fcst[["m2x"]][,c("fcst","lower","upper")]

data <- data.frame(
  quarter = as.Date("2016-02-15") + (0:22)*90,
  value = as.vector(end.vars[,"m2x"])
)
data[c("forecast","lower","upper")] <- NA
data[(dim(data)[1]-3):dim(data)[1],c("forecast","lower","upper")] <- pred3
mape3 <- round(mape(test,pred3[,1]),3)

ggplot(data, aes(x=quarter)) +
  geom_line(aes(y=value)) +
  geom_line(aes(y=upper), color="red", linetype="dashed", size=1) +
  geom_line(aes(y=forecast), color="blue", linetype="longdash", size=1) +
  geom_line(aes(y=lower), color="red", linetype="dashed", size=1) +
  theme_light() + labs(x="", y="") + 
  ggtitle(str_glue("Широкая денежная масса M2X (MAPE = {mape3})"))

cotton <- ex.vars[,"cotton"]
plot(cotton)
acf(cotton)
fit <- auto.arima(cotton, ic="bic")
fit_pred <- forecast(fit,h=13,level=0.95)
plot(forecast(fit,h=13,level=0.95), shaded=FALSE,
     fcol="blue", pi.col="red", main="Вневыборочный прогноз цен на хлопок")

money <- c("m0","m1","m2","m2x")
mains <- c("Наличная денежная масса М0","Узкая денежная масса М1",
           "Широкая денежная масса в нац. определении М2","Широкая денежная масса М2Х")
q = as.Date("2016-02-15") + (0:35)*92
for (i in 1:4) {
  m <- money[i]
  
  train <- end.vars[,c("gdp_real","cpi","lending_rate",m)]
  train_x <- window(ex.vars[,"cotton"], start=c(2016,1), end=c(2021,3))
  train_x <- matrix(train_x,length(train_x),1)
  test_x <- matrix(fit_pred$mean,13,1)
  
  mod <- VAR(train, p=1, type="both", season=4, exogen=train_x)
  pred <- predict(mod, n.ahead=13, dumvar=test_x)$fcst[[m]][,c("fcst","lower","upper")]
  write.csv(pred,paste0(m,".csv"))
  data <- data.frame(
    quarter = q,
    value = c(as.vector(end.vars[,m]),rep(NA,13))
  )
  data[c("forecast","lower","upper")] <- NA
  data[(dim(data)[1]-12):dim(data)[1],c("forecast","lower","upper")] <- pred
  plt <- ggplot(data, aes(x=quarter)) +
    geom_line(aes(y=value)) +
    geom_line(aes(y=upper), color="red", linetype="dashed", size=1) +
    geom_line(aes(y=forecast), color="blue", linetype="longdash", size=1) +
    geom_line(aes(y=lower), color="red", linetype="dashed", size=1) +
    theme_light() + labs(x="", y="") + 
    ggtitle(mains[i])
  print(plt)
}
