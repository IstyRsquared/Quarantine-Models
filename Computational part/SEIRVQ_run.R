## Quarantine model 
## Name: Isty Rysava
## Date: 08/02/2023
## Code: an example script to SEIR.tauleap model (NB) to test intervention scenarios for R0=1

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

### Libraries 
source("R/TauLeap_SEIRfc_vacc_explicit_updated.R")
source("R/quarantine_params.R")

### Parameters and initial conditions 
## Simulation
end.time <- 52*10
nsim <- 1000

## Disease
R0s <- seq(1, 2, 0.1) 
R0 <- R0s[1]
sqcs <- 1:3 
vcs <-  seq(0, 1, 0.05)
vcs <- ifelse(is.infinite(-log(1-vcs)/53), 3/53, -log(1-vcs)/53)

params_grid <- expand.grid(list(R0 = R0, # reproductive number
                                sqc = sqcs, # quarantine scenarios
                                vc = vcs)) # vaccination scenarios
# Pops
N <- 200000
I <- round(49/56/0.05, digits=0) # assuming surveillance detected 5% of all cases in 13 months=56 weeks
E <- round(4*I, digits=0)

allout <- vector("list", length=nrow(params_grid))

for(idx in 1:nrow(params_grid)){
  ### params
  parameters <- params_list[[params_grid$sqc[idx]]] 
  parameters["R0"] <-  params_grid$R0[idx]
  parameters["Rprob"] <- parameters["size"]/(parameters["R0"]+parameters["size"])
  parameters["vc"] <-  params_grid$vc[idx]
  
  ### inits
  V <- ifelse(parameters["vc"]==0, 0, round((-log(parameters["vc"]*53)-1)*N, digits=0))
  V=0.25*N
  V.comp <- rowSums(rmultinom(n=V, size=1, prob=c(50, 30, 20)))
  V1 <- V.comp[1]; V2 <- V.comp[2]; V3 <- V.comp[3]
  S <- N - I - E - V
  
  initials <- c(S=as.numeric(S), V1=V1, V2=V2, V3=V3, E=E, I=I, Qs=0, Qer=0, Qeb=0, Qi=0,
                Rill=0, Sh=1324935, Eh=0, Ih=0,Rh=0, Vhs=0, Vhe=0) 
  
  ### run the model
  deadD <- matrix(NA, ncol=nsim, nrow=end.time+1)
  expD <- matrix(NA, ncol=nsim, nrow=end.time+1) 
  infD <- matrix(NA, ncol=nsim, nrow=end.time+1) 
  deadH <- matrix(NA, ncol=nsim, nrow=end.time+1)
  
  for(sim in 1:nsim){
    sierTL.out <- SEIR.tauleap(init = initials, pars = parameters,
                               end.time = end.time, tau=1) # 1/52 weekly steps
    pops.temp <- as.matrix(sierTL.out$counts)
    res.temp <- as.matrix(sierTL.out$results)
    colnames(pops.temp) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh",
                             "Ih", "Rh", "Vhs", "Vhe")
    deadD[,sim] <- pops.temp[,"Rill"]
    expD[,sim] <- res.temp[,"E"]
    infD[,sim] <- res.temp[,"I"]
    deadH[,sim] <- pops.temp[,"Rh"]
    rm( sierTL.out)
  }
  
  ### extract info & save
  allout[[idx]] <- list(deadD=deadD, expD=expD, infD=infD, deadH=deadH)
  print(idx)
  print(params_grid[idx,])
  
}

# saveRDS(allout, paste0("output/MS_sim_runs_R0", R0, ".Rdata")) 

