## Quarantine model 
## Name: Isty Rysava
## Date: 07/11/2021
## Code: Theoretical part - SEI(R)VQ model & stability analysis (alternative route)
## varying R0 and q 

rm(list=ls())
setwd("~/Dropbox/HongKong")

### Libraries 
require(deSolve)
require(rootSolve)
library(harrypotter)
library("grDevices")
source("R/thesis_params.R")
source("R/thesis_vaccination_params.R")
source("R/detSEIRVQ_rabies.R")

## Initial conditions
I <- I
E <- E
R <- 0
Q <- 0
V <- V1+V2+V3
(N <- S+E+I+V+Q)

## Parameters (by week)
b <- as.numeric(params1["b"]) # birth rate
d <- as.numeric(params1["d"]) # death rate
sigma <- as.numeric(params1["sigma"]) # incubation rate
gamma <- as.numeric(params1["gamma"]) # recovery rate
vc <- -log(1-0.3)/52 # proportion of susceptible dogs vaccinated each week = 30% a year, converted to a rate
wn <- 1/(52*3) # immunity loss each week; vaccine lonegevity ~ 3 yeras
q <- -log(1-0.05)/1 # proportion of infected dogs quarantined each week 5%, converted to a rate
del <- (1/2)*7/7 # delay before quarantined
tau <- gamma/(1-(gamma*del)) # removal rate from the quarantined class
K <- 500000 # carrying capacity
delta <- (b-d)/K
# R0 <- 1.5 # R0
# beta <- ((d+gamma+q+delta*N)*(d+sigma+delta*N)*R0)/sigma # transmisison rate 

##---------------------------------STABILITY: Endemic Eq ---------------------------------##
## Loop through a range of R0 and qs values
R0s <- seq(1,2,.05) # 21
qs <- seq(0, 0.95, 0.05) # 17
N

eigens.all <- c()
periods.all <- c()

Sstar.all <- c()
Estar.all <- c()
Istar.all <- c()
Qstar.all <- c()
Vstar.all <- c()

for(R0 in R0s){
  beta <- ((d+gamma+delta*N)*(d+sigma+delta*N)*R0)/sigma
  # beta <- ((d+gamma+delta)*(d+sigma+delta)*R0)/sigma
  
  for(q.temp in qs){
    q <- -log(1-q.temp)/1
    
    print(paste0("R0:", R0))
    print(paste0("beta:", beta))
    print(paste0("q:", q))
    
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
    # dR <- expression(gamma*I + tau*Q + d*S + d*E + d*I + d*V + d*Q + delta*N*S + delta*N*E + delta*N*I + delta*N*V + delta*N*Q)
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
    eigens.all <- c(eigens.all, eigens.tmp) # note: complex numbers always come in conjugate pairs
    
  }
}


## Get eigenvalues 
# eigenvalues 
eigens.final <- as.data.frame(matrix(eigens.all, ncol=5, byrow=T))
colnames(eigens.final) <- c("eigS", "eigE", "eigI", "eigQ", "eigV")
sign <- rep(1, nrow(eigens.final)) # will be 1 for at least one POS, and 0 for all NEG
for(i in 1:nrow(eigens.final)){
  if(all(eigens.final[i,]<0)){
    sign[i] <- 0 
  }
}
eigens.final$sign <- sign
eigens.final$R0 <- rep(R0s, each=length(qs))
eigens.final$qP <- rep(qs, times=length(R0s))
eigens.final$q <- -log(1-eigens.final$qP)/1

eigens.final$Sstar <- Sstar.all
eigens.final$Estar <- Estar.all
eigens.final$Istar <- Istar.all
eigens.final$Qstar <- Qstar.all
eigens.final$Vstar <- Vstar.all

head(eigens.final)
which(eigens.final$sign==0) # STABLE (note: SEIRVQ will never be stabe bc eigQ is alwas 0)
eigens.final[which(eigens.final$sign==0),]
which(eigens.final$sign==1) # UNSTABLE
eigens.final[which(eigens.final$sign==1),]

## Get periods
eigens.final$periods_wks <- periods.all
head(eigens.final); tail(eigens.final)
eigens.final[which(eigens.final$periods_wks!=Inf),] 
sort(as.numeric(as.character(eigens.final$periods_wks[which(eigens.final$periods_wks!=Inf)]))) 
# write.csv(eigens.final, "~/Dropbox/HongKong/output/SEIRVQrabies_Mike.csv", row.names=F)

## Plot isoclines for each eigenvalue and a heatmap for the sign
pal <- hp(n = 10, house = "Slytherin")
plot(1:10, 1:10, col=pal, pch=16, cex=3)
colgreen <- pal[5]
collightgreen <- pal[2]

pal <- hp(n=8, option = "LunaLovegood")
plot(1:8, 1:8, col=pal, pch=16, cex=3)
colpink <- pal[6]

col <- colorRampPalette(c(colgreen, colpink))(30)
col.more <- colorRampPalette(c(collightgreen, colpink))(60)
plot(1:30, 1:30, col=col, pch=16, cex=3)
plot(1:60, 1:60, col=col.more, pch=16, cex=3)

# plot eigenvalues
pdf("~/Dropbox/HongKong/figs/SEIRVQrabies_Sign_Mike.pdf", width=9, height=5) 
par(mfrow=c(1, 1))
z=matrix(eigens.final$sign, ncol=length(y), nrow=length(x))
image(x=x, y=y, z=z, xlab="q", ylab="R0", col=col, main="Stability", cex.axis=.6, cex.main=.75, cex.lab=.7)
dev.off()

pdf("~/Dropbox/HongKong/figs/SEIRVQrabies_Periods_wks_Mike.pdf", width=6, height=5)
par(mfrow=c(1, 1))
eigens.final$periods_wks[which(eigens.final$periods_wks==Inf)] <- 0
eigens.final$periods_wks <- as.numeric(as.character(eigens.final$periods_wks))
x=eigens.final$qP[1:length(qs)]
y=unique(eigens.final$R0)
z=matrix(eigens.final$periods_wks, ncol=length(y), nrow=length(x))
image(x=x, y=y, z=z, xlab="qP", ylab="R0", col=col.more, main="Periodicity in weeks", cex.axis=.6, cex.main=.75, cex.lab=.7)
contour(x=x, y=y, z=z, add=T, lwd=.2, labcex=.5)
dev.off()


##---------------------------------STABILITY: Disease-free Eq ---------------------------------##
I <- 0
E <- 0
R <- 0
Q <- 0
V <- V1+V2+V3
(N <- S+E+I+V+Q)

## Loop through a range of R0 and qs values
R0s <- seq(1,2,.05) # 21
qs <- seq(0, 0.95, 0.05) # 17
N

eigens.all <- c()
periods.all <- c()

Sstar.all <- c()
Estar.all <- c()
Istar.all <- c()
Qstar.all <- c()
Vstar.all <- c()

for(R0 in R0s){
  beta <- ((d+gamma+delta*N)*(d+sigma+delta*N)*R0)/sigma
  # beta <- ((d+gamma+delta)*(d+sigma+delta)*R0)/sigma
  
  for(q.temp in qs){
    q <- -log(1-q.temp)/1
    
    print(paste0("R0:", R0))
    print(paste0("beta:", beta))
    print(paste0("q:", q))
    
    ## Populations at Endemic Equilibria
    Sstar <- (b*N*(d+delta*N+wn))/((d+delta*N+vc)*(d+delta*N+wn)-wn*vc)
    Estar <- 0
    Istar <- 0
    Vstar <- (vc*b*N)/((d+delta*N+vc)*(d+delta*N+wn)-wn*vc)
    Qstar <- 0
    
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
    # dR <- expression(gamma*I + tau*Q + d*S + d*E + d*I + d*V + d*Q + delta*N*S + delta*N*E + delta*N*I + delta*N*V + delta*N*Q)
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
    eigens.all <- c(eigens.all, eigens.tmp) # note: complex numbers always come in conjugate pairs
    
  }
}


## Get eigenvalues 
# eigenvalues 
eigens.final <- as.data.frame(matrix(eigens.all, ncol=5, byrow=T))
colnames(eigens.final) <- c("eigS", "eigE", "eigI", "eigQ", "eigV")
sign <- rep(1, nrow(eigens.final)) # will be 1 for at least one POS, and 0 for all NEG
for(i in 1:nrow(eigens.final)){
  if(all(eigens.final[i,]<0)){
    sign[i] <- 0 
  }
}
eigens.final$sign <- sign
eigens.final$R0 <- rep(R0s, each=length(qs))
eigens.final$qP <- rep(qs, times=length(R0s))
eigens.final$q <- -log(1-eigens.final$qP)/1

eigens.final$Sstar <- Sstar.all
eigens.final$Estar <- Estar.all
eigens.final$Istar <- Istar.all
eigens.final$Qstar <- Qstar.all
eigens.final$Vstar <- Vstar.all

head(eigens.final)
which(eigens.final$sign==0) # STABLE (note: SEIRVQ will never be stabe bc eigQ is alwas 0)
eigens.final[which(eigens.final$sign==0),]
which(eigens.final$sign==1) # UNSTABLE
eigens.final[which(eigens.final$sign==1),]

## Get periods
eigens.final$periods_wks <- periods.all
head(eigens.final); tail(eigens.final)
eigens.final[which(eigens.final$periods_wks!=Inf),] 
sort(as.numeric(as.character(eigens.final$periods_wks[which(eigens.final$periods_wks!=Inf)]))) 
# write.csv(eigens.final, "~/Dropbox/HongKong/output/SEIRVQrabies_DiseaseFreeEq.csv", row.names=F)

## Plot isoclines for each eigenvalue and a heatmap for the sign
pal <- hp(n = 10, house = "Slytherin")
plot(1:10, 1:10, col=pal, pch=16, cex=3)
colgreen <- pal[5]
collightgreen <- pal[2]

pal <- hp(n=8, option = "LunaLovegood")
plot(1:8, 1:8, col=pal, pch=16, cex=3)
colpink <- pal[6]

col <- colorRampPalette(c(colgreen, colpink))(30)
col.more <- colorRampPalette(c(collightgreen, colpink))(60)
plot(1:30, 1:30, col=col, pch=16, cex=3)
plot(1:60, 1:60, col=col.more, pch=16, cex=3)

# plot eigenvalues
pdf("~/Dropbox/HongKong/figs/SEIRVQrabies_Sign_DiseaseFreeEq.pdf", width=9, height=5) 
x=eigens.final$qP[1:length(qs)] 
y=unique(eigens.final$R0) 
par(mfrow=c(1, 1))
z=matrix(eigens.final$sign, ncol=length(y), nrow=length(x))
image(x=x, y=y, z=z, xlab="q", ylab="R0", col=col, main="Stability", cex.axis=.6, cex.main=.75, cex.lab=.7)
dev.off()




