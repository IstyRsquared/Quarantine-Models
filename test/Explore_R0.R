## Quarantine model 
## Name: Isty Rysava
## Date: 20/12/22
## Code: explore formualtions for generating offspring cases

rm(list=ls())
# setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

### Init
probrab <- 1 # 1 or 0.49
R0s <- seq(1, 2, 0.1) 
sizes <- seq(1, 2, 0.1)
probs <- seq(0.05, 0.95, 0.05)

### Test
## Poiss
pdf(paste0("figs/R0_explore/Poiss_probrab", probrab, ".pdf"), width=4, height=8)
par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
for(R0 in R0s){
  offspring <- rpois(n=1000, lambda=R0)
  ## add prob rab
  hist(offspring, breaks=seq(0, max(offspring), 1), main=paste0("R0=", R0),
       cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
}
dev.off()

## Neg Binom: through mean
pdf(paste0("figs/R0_explore/NegBin_mean_probrab", probrab, ".pdf"), width=4, height=8)
par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
for(R0 in R0s){
  for(size.temp in sizes){
    prob.temp <- size.temp/(R0+size.temp)
    offspring <- rnbinom(n=1000, size=size.temp, prob=prob.temp)
    ## add prob rab
    hist(offspring, breaks=seq(0, max(offspring), 1), 
         main=paste0("R0=", R0, " size=", size.temp, " prob=", round(prob.temp, digits=2)),
         cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
  }
}
dev.off()

## Neg Binom: size and prob
pdf(paste0("figs/R0_explore/NegBin_size_probrab", probrab, ".pdf"), width=4, height=8)
par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
for(R0 in R0s){
  for(prob.temp in probs){
    offspring <- rnbinom(n=1000, size=R0, prob=prob.temp)
    ## add prob rab
    hist(offspring, breaks=seq(0, max(offspring), 1),
         main=paste0("R0=", R0, " prob=", prob.temp),
         cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
  }
}
dev.off()

############ Draw plots by R0
for(R0 in R0s){
  pdf(paste0("figs/R0_explore/Explore_R0", R0, ".pdf"), width=4, height=8)
  par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
  
  ## 1) Poiss
  offspring <- rpois(n=1000, lambda=R0)
  hist(offspring, breaks=seq(0, max(offspring), 1), main=paste0("type=Poiss"),
       cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
  
  ## 2) Neg Binom: through mean
  for(size.temp in sizes){
    prob.temp <- size.temp/(R0+size.temp)
    offspring <- rnbinom(n=1000, size=size.temp, prob=prob.temp)
    ## add prob rab
    hist(offspring, breaks=seq(0, max(offspring), 1), 
         main=paste0("type=NBmean", " size=", size.temp, " prob=", round(prob.temp, digits=2)),
         cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
  }
  
  ## 3) Neg Binom: size and prob
  for(prob.temp in probs){
    offspring <- rnbinom(n=1000, size=R0, prob=prob.temp)
    ## add prob rab
    hist(offspring, breaks=seq(0, max(offspring), 1),
         main=paste0("type=NBsize", " prob=", prob.temp),
         cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
  }
  dev.off()
}