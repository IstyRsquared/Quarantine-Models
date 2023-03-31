## Quarantine model 
## Name: Isty Rysava
## Date: 31/03/2023
## Code: Statistical analysis of vacination and quarantine effect on the number of cases

rm(list=ls())

# setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

# library(harrypotter)
library(tidyverse)
library(MASS)
library(broom)

## Read in data and intit
final_frame_box <- read.csv("output/MS_monthly_infection_boxplot.csv")
R0s <- seq(1, 2, 0.1)

#################################################################################################################################################
### Get stats
# NOTE: check if POiss would be a better fit

statD_res <- matrix(NA, nrow=length(R0s), ncol=7)
statH_res <- matrix(NA, nrow=length(R0s), ncol=7)
statD_res[,1] <- statH_res[,1] <- R0s

for(idx in 1:length(R0s)){
  # subset & get total infection
  Rsubset <- filter(final_frame_box, R0==R0s[idx])
  Rsubset$allinfD <- Rsubset$expD + Rsubset$infD 
  
  # run model for dogs & get sig and ests
  mod1 <- glm.nb(allinfD ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
  statD_res[idx,2:4] <- as.numeric(tidy(mod1)$p.value[2:4])
  statD_res[idx,5:7] <- as.numeric(exp(tidy(mod1)$estimate)[2:4])
  # the incident rate for Q2 is 0.85 times the incident rate for the reference group (Q1)
  # the incident rate for Q3 is 0.67 times the incident rate for the reference group 
  # the percent change in the incident rate of cases is a 0.35% decrease for every unit increase in vacc

  # run  mod for humans & get sig and ests
  mod2 <- glm.nb(deadH ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
  statH_res[idx,2:4] <- as.numeric(tidy(mod2)$p.value[2:4])
  statH_res[idx,5:7] <- as.numeric(exp(tidy(mod2)$estimate)[2:4])
}

colnames(statD_res) <- colnames(statH_res) <- c("R0", "pvalVacc", "pvalQ2", "pvalQ3", "VaccQ3", "estQ2", "estQ3")
colnames(statD_res) <- data.frame(statD_res)
colnames(statH_res) <- data.frame(statH_res)

# write.csv(statD_res, "output/statD_res.csv", row.names=F)
# write.csv(statH_res, "output/statH_res.csv", row.names=F)

#################################################################################################################################################
### Plot predicted
Rsubset <- filter(final_frame_box, R0==1.3)
Rsubset$allinfD <- Rsubset$expD + Rsubset$infD 

# run model for dogs
mod1 <- glm.nb(allinfD ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
summary(mod1)
newdata <- expand.grid(Vaccination = c(0, 0.25, 0.5, 0.75), Quarantine = factor(1:3))
newdata <- cbind(newdata, predict(mod1, newdata, type = "link", se.fit=TRUE))
newdata <- within(newdata, {
  Cases <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(Vaccination, Cases)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = Quarantine), alpha = .25) +
  geom_line(aes(colour = Quarantine), size = 0.1) +
  labs(x = "Vaccination", y = "Monthly no. of infected dogs (E+I)")

# run  mod for humans 
mod2 <- glm.nb(deadH ~ Vaccination + as.factor(Quarantine), Rsubset, na.action=na.exclude, maxit=1000, link=log)
summary(mod2)
newdata <- expand.grid(Vaccination = c(0, 0.25, 0.5, 0.75), Quarantine = factor(1:3))
newdata <- cbind(newdata, predict(mod2, newdata, type = "link", se.fit=TRUE))
newdata <- within(newdata, {
  Cases <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

ggplot(newdata, aes(Vaccination, Cases)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = Quarantine), alpha = .25) +
  geom_line(aes(colour = Quarantine), size = 0.5) +
  labs(x = "Vaccination", y = "Monthly no. of human fatalities due to rabies")


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
