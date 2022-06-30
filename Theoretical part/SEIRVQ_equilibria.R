## Quarantine model 
## Name: Isty Rysava
## Date: 06/29/2022
## Code: Theoretical part - SEI(R)VQ model & stability analysis
## varying R0, q, and vc 

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

## Libraries 
require(deSolve)
require(rootSolve)
library(harrypotter)
library("grDevices")

source("R/model_params.R")

## Parameters (by week)
b <- as.numeric(params1["b"]) # birth rate
d <- as.numeric(params1["d"]) # death rate
sigma <- as.numeric(params1["sigma"]) # incubation rate
gamma <- as.numeric(params1["gamma"]) # recovery rate
wn <- 1/(52*3) # immunity loss each week; vaccine lonegevity ~ 3 yeras
del <- (1/2)*7/7 # delay before quarantined
tau <- gamma/(1-(gamma*del)) # removal rate from the quarantined class
K <- 500000 # carrying capacity
delta <- (b-d)/K

##---------------------------------STABILITY: Endemic Eq ---------------------------------##
## Set up parametr grid to explore
R0s <- seq(1, 2, 0.1) 
qs <- seq(0, 1, 0.05) 
qs <- ifelse(is.infinite(-log(1-qs)/1), qs, -log(1-qs)/1)
vcs <-  seq(0, 1, 0.05)
vcs <- ifelse(is.infinite(-log(1-vcs)/1), vcs, -log(1-vcs)/1)

params_grid <- expand.grid(list(R0 = R0s, # reproductive number
                                q = qs, # quarantine rate
                                vc = vcs)) # vaccination rate

## Initialize
N <- 200000

eigens.all <- matrix(NA, ncol=5, nrow=nrow(params_grid))
periods.all <- c()

Sstar.all <- c()
Estar.all <- c()
Istar.all <- c()
Qstar.all <- c()
Vstar.all <- c()

## Loop through parametr space
for(idx in 1:nrow(params_grid)){
  R0 <- params_grid$R0[idx]
  q <- params_grid$q[idx]
  vc <- params_grid$vc[idx]
  beta <- ((d+gamma+delta*N)*(d+sigma+delta*N)*R0)/sigma
  
  ## Populations at Endemic Equilibria
  Sstar <- (N*(d+delta*N+sigma)*(d+delta*N+gamma+q))/(beta*sigma)
  Estar <- ((b*N)/(d+delta*N+sigma))  - (N*(d+delta*N+vc)*(d+delta*N+gamma+q))/(beta*sigma) + (wn*N*vc*(d+delta*N+gamma+q))/(beta*sigma*(d+delta*N+wn))
  Istar <- (b*N*sigma)/((d+delta*N+sigma)*(d+delta*N+gamma+q)) - (N*(d+delta*N+vc))/beta + (wn*N*vc)/(beta*(d+delta*N+wn))
  Vstar <- (vc*N*(d+delta*N+sigma)*(d+delta*N+gamma+q))/((beta*sigma)*(d+delta*N+wn))
  Qstar <- q/(d+delta*N+tau)*((b*N*sigma)/((d+delta*N+sigma)*(d+delta*N+gamma+q)) - (N*(d+delta*N+vc))/beta + (wn*N*vc)/(beta*(d+delta*N+wn)))
  
  eq1 <- list(S=Sstar, E=Estar, I=Istar, Q=Qstar, V=Vstar)
  Sstar.all <- c(Sstar.all, Sstar)
  Estar.all <- c(Estar.all, Estar)
  Istar.all <- c(Istar.all, Istar)
  Qstar.all <- c(Qstar.all, Qstar)
  Vstar.all <- c(Vstar.all, Vstar)
  
  ## Equations
  dS <- expression(b*N - d*S - delta*N*S - beta*S*I/N - vc*S + wn*V)
  dE <- expression(beta*S*I/N - d*E - delta*N*E - sigma*E)
  dI <- expression(sigma*E - d*I - delta*N*I - gamma*I - q*I)
  dQ <- expression(q*I - d*Q - delta*N*Q - tau*Q)
  dV <- expression(vc*S - d*V - delta*N*V - wn*V)
  
  j11 <- D(dS, "S"); j12 <- D(dS, "E"); j13 <- D(dS, "I"); j14 <- D(dS, "Q"); j15 <- D(dS, "V")
  j21 <- D(dE, "S"); j22 <- D(dE, "E"); j23 <- D(dE, "I"); j24 <- D(dE, "Q"); j25 <- D(dE, "V")
  j31 <- D(dI, "S"); j32 <- D(dI, "E"); j33 <- D(dI, "I"); j34 <- D(dI, "Q"); j35 <- D(dI, "V")
  j41 <- D(dQ, "S"); j42 <- D(dQ, "E"); j43 <- D(dQ, "I"); j44 <- D(dQ, "Q"); j45 <- D(dQ, "V")
  j51 <- D(dV, "S"); j52 <- D(dV, "E"); j53 <- D(dV, "I"); j54 <- D(dV, "Q"); j55 <- D(dV, "V")
  
  ## Evaluate Jacobian at equilibrium
  J <- with(data=eq1, expr = matrix(c(eval(j11), eval(j12), eval(j13), eval(j14), eval(j15), 
                                      eval(j21), eval(j22), eval(j23), eval(j24), eval(j25), 
                                      eval(j31), eval(j32), eval(j33), eval(j34), eval(j35), 
                                      eval(j41), eval(j42), eval(j43), eval(j44), eval(j45), 
                                      eval(j51), eval(j52), eval(j53), eval(j54), eval(j55)), 
                                    nrow=5, byrow=T))
  
  ## Calculate eigenvalues 
  eigens.tmp <- as.numeric(as.character(Re(eigen(J)$values)))
  de.ind <- which.max(Re(eigen(J)$values))
  periods.all <- c(periods.all, as.character(2*pi/Im(eigen(J, only.values=T)$values[de.ind])))
  eigens.all[idx,] <- eigens.tmp # note: complex numbers always come in conjugate pairs
  
}

## Summarize 
# eigenvalues sign
head(eigens.all)
sign <- rep(1, nrow(eigens.all)) # will be 1 for at least one POS, and 0 for all NEG
for(i in 1:nrow(eigens.all)){
  if(all(eigens.all[i,]<0)){
    sign[i] <- 0 
  }
}

# infection
infection <- Istar.all + Estar.all 

# combine in df
params <- expand.grid(list(R0 = seq(1, 2, 0.1), # reproductive number
                           q = seq(0, 1, 0.05) , # quarantine rate
                           vc = seq(0, 1, 0.05))) # vaccination rate

sign_df <- cbind(params, variable=rep("stability", length(sign)), value=as.numeric(sign))
period_df <- cbind(params, variable=rep("period", length(periods.all)), value=as.numeric(periods.all))
inf_df <- cbind(params, variable=rep("infection", length(infection)), value=as.numeric(infection))

final_df <- rbind(sign_df, period_df, inf_df)
head(final_df); tail(final_df)
# write.csv(final_df, "output/SEIRVQrabies_EndemicEquilibrium.csv", row.names=F)



