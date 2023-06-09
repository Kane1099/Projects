#Ефимов Александр э541

library("dplyr")
library("psych")
library("lmtest")
library("glmnet")
library("ggplot2")
library("gplots")
library("car")
library("zeallot")
library("hdm")
library("pbapply")
library("grf")
library("caret")

set.seed(1234)

data<-read.csv("data2021_HA3upd.csv", header=TRUE, dec=".", sep=";")
data$price_std<-as.numeric(data$price_std)
data$final_rating_mean<-as.numeric(data$final_rating_mean)
data<-na.omit(data)

# Проверка наличия мультиколлинеарности
m0 <- lm(total_profit ~ ., data = data[-c(1,39)])
summary(m0)
m0_vif <- vif(m0)
m0_vif <- data.frame(val=m0_vif)

m0_vif %>%
  ggplot(aes(reorder(rownames(m0_vif),val),val)) +
  geom_col(fill="#f68060", alpha=.6, width=.4) +
  geom_text(aes(label = round(val,2)), hjust = -0.05) +
  geom_hline(yintercept = 10, linetype = 'dashed') +
  coord_flip() +
  xlab("") + ylab("VIF") +
  theme_bw()
  
# Double lasso
X0 <- model.matrix(data = data[-c(1,2,39)], total_profit ~ .)
Y0 <- data$total_profit
X1 <- model.matrix(data = data[-c(1,2,40)], buy ~ .)
Y1 <- data$buy
W <- data$T
Eff0 <- rlassoEffect(X0, Y0, W, method = "double selection")
coefficients(Eff0)
Eff1 <- rlassologitEffect(X1, Y1, W)
coefficients(Eff1)

results <- data.frame(tau=c(coefficients(Eff0),coefficients(Eff1)),
                      lower=c(confint(Eff0)[1],confint(Eff1)[1]),
                      upper=c(confint(Eff0)[2],confint(Eff1)[2]))
rownames(results) <- c("total_profit","buy")

ggplot(results, aes(x = rownames(results), y = tau, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_pointrange() +
  geom_text(label = round(results$tau,3), vjust = -0.75) +
  labs(x = "", y = "Treatment effect") +
  coord_flip()

# Делим выборку на обучение и тест (50/50)
train.index <- createDataPartition(data$T, p = .5, list = FALSE)
train <- data[ train.index,]
test  <- data[-train.index,]

# Причинный лес
tau.forest <- causal_forest(train[c(4,14,15,19:30)], train$total_profit, train$T, num.trees = 500, seed = 1234)
tau.hat <- predict(tau.forest, test[c(4,14,15,19:30)], estimate.variance = TRUE)

sigma.hat <- sqrt(tau.hat$variance.estimates)
tau_hat <- tau.hat$predictions
CI <- list(L = tau.hat$predictions - qnorm(0.95) * sigma.hat,
          U = tau.hat$predictions + qnorm(0.95) * sigma.hat)

D <- data.frame(tau_pr=tau_hat, l_pr=CI$L, u_pr=CI$U)

tau.forest1 <- causal_forest(train[c(4,14,15,19:30)], train$buy, train$T, num.trees = 500, seed = 1234)
tau.hat1 <- predict(tau.forest1, test[c(4,14,15,19:30)], estimate.variance = TRUE)

sigma.hat1 <- sqrt(tau.hat1$variance.estimates)
tau_hat1 <- tau.hat1$predictions
CI1 <- list(L = tau.hat1$predictions - qnorm(0.95) * sigma.hat1,
           U = tau.hat1$predictions + qnorm(0.95) * sigma.hat1)

D$tau_buy <- tau_hat1
D$l_buy <- CI1$L
D$u_buy <- CI1$U

par(mfrow=c(1,2))
hist(D$tau_pr, breaks=30, col=rgb(1,0,0,0.5) , xlab="total profit" , ylab="" , main="" )
hist(D$tau_buy, breaks=30, col=rgb(0,0,1,0.5) , xlab="buy" , ylab="" , main="")

# Решаем, какому наблюдению какую наценку назначить
# Если эффект значим и нижняя граница > 0, то выбираем T = 1
attach(D)
D11 <- sum(l_pr <= 0 & l_buy <= 0)
D01 <- sum(l_pr > 0 & l_buy <= 0)
D10 <- sum(l_pr <= 0 & l_buy > 0)
D00 <- sum(l_pr > 0 & l_buy > 0)

final <- data.frame(c(D00,D10),c(D01,D11))
colnames(final) <- c("buy effect > 0","buy effect = 0")
rownames(final) <- c("total profit effect > 0","total profit effect = 0")
final

# Либо можно просто смотреть на значение эффекта без учета значимости
D11 <- sum(tau_pr <= 0 & tau_buy <= 0)
D01 <- sum(tau_pr > 0 & tau_buy <= 0)
D10 <- sum(tau_pr <= 0 & tau_buy > 0)
D00 <- sum(tau_pr > 0 & tau_buy > 0)

final <- data.frame(c(D00,D10),c(D01,D11))
colnames(final) <- c("buy effect > 0","buy effect < 0")
rownames(final) <- c("total profit effect > 0","total profit effect < 0")
final