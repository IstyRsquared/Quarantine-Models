## Quarantine model 
## Name: Isty Rysava
## Date: 31/03/2023
## Code: Statistical analysis of vaccination and quarantine effect on the number of cases

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

# library(harrypotter)
library(tidyverse)
library(MASS)
library(broom)

## Read in data and intit
final_frame_box <- read.csv("output/MS_monthly_infection_boxplot_burnout_10yrs.csv")
final_frame_box$R0 <- as.character(final_frame_box$R0)
R0s <- seq(1, 2, 0.1)

#################################################################################################################################################
### Get the statistics
statD_res <- matrix(NA, nrow=length(R0s), ncol=8)
statH_res <- matrix(NA, nrow=length(R0s), ncol=8)
statD_res[,1] <- statH_res[,1] <- R0s

for(idx in 1:length(R0s)){
  # subset & get total infection
  Rsubset <- subset(final_frame_box, R0==as.character(R0s[idx]))
  Rsubset$allinfD <- Rsubset$expD + Rsubset$infD 
  
  # run model for dogs & get sig and ests
  mod1 <- glm.nb(allinfD ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
  statD_res[idx,2:4] <- as.numeric(tidy(mod1)$p.value[2:4])
  statD_res[idx,5:8] <- as.numeric(exp(tidy(mod1)$estimate)[1:4]) 
  ## for R0=1.3
  # intercept 15.8976936
  # the incident rate for Q2 is 0.9144578 times the incident rate for the reference group (Q1)
  # i.e., 15.8976936*0.9144578=14.53777
  # the incident rate for Q3 is 0.8702091 times the incident rate for the reference group (Q1)
  # i.e., 15.8976936*0.8702091=13.83432
  # the percent change in the incident rate of cases is a 0.43% decrease for every unit increase in vacc
  # i.e., 15.8976936*(0.4384563^0.5) [intercept*(coeff^vaccunits)]
      
  # run  mod for humans & get sig and ests
  mod2 <- glm.nb(deadH ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
  statH_res[idx,2:4] <- as.numeric(tidy(mod2)$p.value[2:4])
  statH_res[idx,5:8] <- as.numeric(exp(tidy(mod2)$estimate)[1:4])
  
  print(idx)
}

statD_res <- data.frame(statD_res)
statH_res <- data.frame(statH_res)
colnames(statD_res) <- colnames(statH_res) <- c("R0", "pvalVacc", "pvalQ2", "pvalQ3", "intercept", 
                                                "estVacc", "estQ2", "estQ3")

head(statD_res)
write.csv(statD_res, "output/statD_res.csv", row.names=F)
head(statH_res)
write.csv(statH_res, "output/statH_res.csv", row.names=F)

#################################################################################################################################################
# ### Plot predicted
# Rsubset <- filter(final_frame_box, R0=="1.3")
# Rsubset$allinfD <- Rsubset$expD + Rsubset$infD 
# 
# # CALCULATE PREDICTION INTERVAL LIKE IN TS!
# # run model for dogs
# mod1 <- glm.nb(allinfD ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
# summary(mod1)
# newdata <- expand.grid(Vaccination = c(0, 0.25, 0.5, 0.75), Quarantine = factor(1:3))
# newdata <- cbind(newdata, predict(mod1, newdata, type = "link", se.fit=TRUE))
# newdata <- within(newdata, {
#   Cases <- exp(fit)
#   LL <- exp(fit - 1.96 * se.fit)
#   UL <- exp(fit + 1.96 * se.fit)
# })
# 
# ggplot(newdata, aes(Vaccination, Cases)) +
#   geom_ribbon(aes(ymin = LL, ymax = UL, fill = Quarantine), alpha = .25) +
#   geom_line(aes(colour = Quarantine), size = 0.1) +
#   labs(x = "Vaccination", y = "Monthly no. of infected dogs (E+I)")
# 
# # run  mod for humans 
# mod2 <- glm.nb(deadH ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
# summary(mod2)
# newdata <- expand.grid(Vaccination = c(0, 0.25, 0.5, 0.75), Quarantine = factor(1:3))
# newdata <- cbind(newdata, predict(mod2, newdata, type = "link", se.fit=TRUE))
# newdata <- within(newdata, {
#   Cases <- exp(fit)
#   LL <- exp(fit - 1.96 * se.fit)
#   UL <- exp(fit + 1.96 * se.fit)
# })
# 
# ggplot(newdata, aes(Vaccination, Cases)) +
#   geom_ribbon(aes(ymin = LL, ymax = UL, fill = Quarantine), alpha = .25) +
#   geom_line(aes(colour = Quarantine), size = 0.5) +
#   labs(x = "Vaccination", y = "Monthly no. of human fatalities due to rabies")
# 

# ## Calculate human cases
# # low: 4, high: 6
# idx=6
# stats_mthly <- final_list_ts[[idx]]
# dogs <- stats_mthly[[1]]
# humans <- stats_mthly[[4]]
# sum(humans$mean) 
# # R0=1.3, low: 65.076, high: 47.129
# # (65.076 - 47.129) / (65.076/100) # 28%
# # R0=1.2, low: 60.555, high: 45.101
# # (60.555 - 45.101) / (60.555/100) # 26%
# 
# sum(dogs$mean) 
# # R0=1.3, low: 958.651, high: 817.13
# # (958.651 - 817.13) / (958.651/100) # 15%
# # R0=1.2, low: 886.805, high: 776.412
# # (886.805 - 776.412) / (886.805/100) # 13%
