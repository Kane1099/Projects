# Ефимов Александр э-541

library(stargazer)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(plm)
library(lmtest)
library(readxl)
library(regclass)
library(tidyverse)
library(tableone)
library(pglm)
library(MatchIt)
library(Matching)
library(lmtest)
library(sandwich)
library(survival)

data <- read_xlsx("wage_gap.xlsx")

# Panel models
panel_data <- data[c("n","t","cpi_w","SMSA_central","AFQT2","self_conf",
                     "education","years","woman","black","hispanic",
                     "fam_size","married","union","promotion","risk")]

summary(panel_data$cpi_w)
panel_data$cpi_w <- ifelse(panel_data$cpi_w==0,0,log(panel_data$cpi_w))
xlabs <- colnames(panel_data)[-(1:3)]
fmla <- as.formula(paste("cpi_w ~ ", paste(xlabs, collapse = "+")))

pool <- plm(fmla, index = c("n","t"),
            data = panel_data, model = "pooling")

fe_oneway <- plm(fmla, index = c("n","t"),
            data = panel_data, model = "within", effect = "individual")

fe_twoway <- plm(fmla, index = c("n","t"),
                 data = panel_data, model = "within", effect = "twoways")

re <- plm(fmla, index = c("n","t"),
                 data = panel_data, model = "random", effect = "twoways")

waldtest(fe_oneway, fe_twoway)

cse = function(reg) {
  rob = sqrt(diag(vcovHC(reg, type = "HC1")))
  return(rob)
}
clse = function(reg) { 
  G = length(unique(index(reg,"id")))
  N = length(index(reg,"id"))
  dfa = (G/(G - 1))   
  rob = sqrt(diag(dfa*vcovHC(reg, method="arellano", type = "HC1", 
                             cluster = "group")))
  return(rob)
}

pFtest(fe_oneway, pool) # fixed better than pool
plmtest(re) # random better than pool
phtest(fe_oneway, re) # fixed better than random

stargazer(pool, fe_oneway, fe_twoway, re, type = "text", #out = "panel_table.doc",
          column.labels=c("Pooled", "FE-oneway", "FE-twoway", "RE"),
          se = list(cse(pool), clse(fe_oneway), clse(fe_twoway), clse(re)),
          model.numbers = FALSE, header = FALSE,
          dep.var.labels = "log(cpi.w)",
          #covariate.labels = c("SMSA.central", gsub("_", ".", xlabs)),
          omit.stat= "LL", no.space=TRUE, align = TRUE, df = FALSE)

# Matching
match_data <- data[c("n", "t", "cpi_w","SMSA_central", "fam_size", "AFQT2", "education",
                   "HGT_father", "HGT_mother", "self_conf", "size_of_firm", "risk")]
match_data <- na.omit(match_data)
match_data$cpi_w <- ifelse(match_data$cpi_w==0,0,log(match_data$cpi_w))
covariates <- c("fam_size", "AFQT2", "education", "HGT_father",
                "HGT_mother", "self_conf", "size_of_firm", "risk")

tmp <- match_data %>% group_by(n) %>%
  mutate(stayed = mean(SMSA_central)==0, moved = any(diff(SMSA_central)>0),
         when_moved = c(0, diff(SMSA_central)))
# treat1 <- tmp %>% filter(moved==TRUE,t==1994)
# control1 <- tmp %>% filter(stayed==TRUE,t==1994)
treat2 <- tmp %>% filter(when_moved==1)
control2 <- tmp %>% filter(stayed==TRUE)

model_est <- function(treat, control) {
  
  lim <- quantile(control$cpi_w, 0.95)
  control <- control[control$cpi_w < lim,]
  df <- rbind(treat, control)
  table1 <- CreateTableOne(vars=covariates,
                           strata = "SMSA_central", data=df, test=TRUE)
  print("Balance of covariates before matching")
  print(table1)
  
  fmla1 <- SMSA_central~fam_size+AFQT2+education+HGT_father+HGT_mother+self_conf+size_of_firm+risk
  model_1 <- matchit(fmla1, data = df)
  matched <- df[as.logical(model_1$weights),]
  table2 <- CreateTableOne(vars=covariates, 
                           strata = "SMSA_central", data=matched, test=TRUE)
  print("Balance of covariates after matching")
  print(table2)
  ps <- glm(fmla1, data = df, family = binomial(link='logit'))
  stargazer(ps, type = "text")
  
  model_1_t <- lm(cpi_w ~ SMSA_central, data = df, weights = model_1$weights)
  stargazer(model_1_t, type = "text", df = FALSE, digits = 3)
  
  VV <- cbind(df[covariates], ps$fitted.values)
  rr<-Match(Y=df$cpi_w,Tr=df$SMSA_central,X=VV)
  print(summary(rr))

}

# model_est(treat1, control1)
model_est(treat2, control2)
