library(rdrobust)
library(ggplot2)
library(stargazer)
library(tidyverse)

data <- read.csv("outputnew.csv", sep = ";", dec = ",")

#1
Y1 <- data[data$treated==1,]$educ
Y0 <- data[data$treated==0,]$educ
mean_diff <- mean(Y1) - mean(Y0)
mean_diff
lm0 <- lm(educ ~ treated, data = data)
summary(lm0)

#3
data$R <- ifelse(data$treated==1,data$dist2cutoff,-data$dist2cutoff)
ggplot(data = data, aes(x = R, y = educ, color = treated)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(data = filter(data, R <= 0), method = "lm") +
  geom_smooth(data = filter(data, R > 0), method = "lm") +
  geom_vline(xintercept = 0)

#4
lm1 <- lm(educ ~ treated * R, data = data)
summary(lm1)

#6
df6 <- data %>% filter(dist2cutoff <= 3000)
lm2 <- lm(educ ~ treated * R, data = df6)
summary(lm2)

stargazer(lm0, lm1, lm2, type = "text", digit = 3, df = FALSE,
          column.labels = c("Mean diff", "incl. R", "incl. R (dist < 3km)"))

#7,9

rdd_est <- function(cutoff, kernels) {
  
  models <- c()
  
  attach(data)
  for (c in cutoff) {
    for (ker in kernels) {
      rd <- rdrobust(educ, R, kernel = ker, c = c)
      res <- paste(round(rd$Estimate[1],3)," (",round(rd$pv[1],3),")", sep = "")
      models <- c(models, res)
    }
  }
  models <- matrix(models, 3, 5)
  models <- cbind(kernels, models)
  results <- as.data.frame(models)
  colnames(results) <- c("kernel", paste("c=", cut, sep = ""))
  return(results)
}

cut <- c(0,-3000,-1000,1000,3000)
kernels <- c("triangular", "epanechnikov", "uniform")
results <- rdd_est(cut, kernels)
results
