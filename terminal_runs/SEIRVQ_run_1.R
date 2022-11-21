## Quarantine model 
## Name: Isty Rysava
## Date: 08/02/18
## Code: runs  SEIR.tauleap model to test intervention scenarios 

rm(list=ls())
# setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

### Libraries 
source("R/TauLeap_SEIRfc_vacc_explicit_updated.R")
source("R/model_params.R")

### Parameters and initial conditions 
## Simulation
end.time <- 52*5
nsim <- 1000

## Disease
R0s <- seq(1, 2, 0.1) 
R0 <- R0s[1]
sqcs <- 1:3 
vc.temp <-  c(0, 0.25, 0.5, 0.75)
vcs <- ifelse(is.infinite(-log(1-vc.temp)/53), vc.temp, -log(1-vc.temp)/53)

params_grid <- expand.grid(list(R0 = R0, # reproductive number
                                sqc = sqcs, # quarantine scenarios
                                vc = vcs)) # vaccination scenarios
# Pops
N <- 200000
I <- round(49/56/0.05, digits=0) # assuming surveillance detected 5% of all cases in 13 months=56 weeks
E <- round(4*I, digits=0)

# allout <- vector("list", length=nrow(params_grid))

for(idx in 1:nrow(params_grid)){
  ### params
  parameters <- params_list[[params_grid$sqc[idx]]] 
  parameters["R0"] <-  params_grid$R0[idx]
  parameters["Rprob"] <- parameters["size"]/(parameters["R0"]+parameters["size"])
  parameters["vc"] <-  params_grid$vc[idx]
  
  ### inits
  V <- ifelse(parameters["vc"]==0, 0, round((1-exp(-parameters["vc"]*53))*N, digits=0))
  V.comp <- rowSums(rmultinom(n=V, size=1, prob=c(50, 30, 20)))
  V1 <- V.comp[1]; V2 <- V.comp[2]; V3 <- V.comp[3]
  S <- N - I - E - V
  
  initials <- c(S=S, V1=V1, V2=V2, V3=V3, E=E, I=I, Qs=0, Qer=0, Qeb=0, Qi=0,
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
  # allout[[idx]] <- list(deadD=deadD, expD=expD, infD=infD, deadH=deadH)
  allout <- list(deadD=deadD, expD=expD, infD=infD, deadH=deadH)
  print(idx)
  print(params_grid[idx,])
  saveRDS(allout, paste0("output/MS_sim_runs_R0", parameters["R0"], 
                         "_sqc", params_grid$sqc[idx],
                         "_vc", params_grid$vc[idx], ".Rdata")) 
}

# saveRDS(allout, paste0("output/MS_sim_runs_R0", R0, ".Rdata")) 

# ### Run the model: testing
# end.time <- 52*5
# sierTL.out <- SEIR.tauleap(init = initials, pars = parameters,
#                            end.time = end.time, tau=1) # 1/52 weekly steps
# out <- sierTL.out$results
# counts <- sierTL.out$counts
# 
# out <- as.matrix(out)
# counts <- as.data.frame(counts)
# colnames(counts) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh",
#                    "Ih", "Rh", "Vhs", "Vhe")
# head(counts); tail(counts)
# counts$E 
# 
# ### Plot
# # dog population
# plot(seq(0, end.time, by=1), out[,2], type="l", ylim=c(0,max(out[,2], na.rm=T))) # sus
# lines(seq(0, end.time, by=1), rowSums(out[,3:5]), col="orange") # all vacc
# 
# # dog vaccination
# plot(seq(0, end.time, by=1), out[,3],  type="l", col="red", ylim=c(0,max(out[,3]))) # V1
# lines(seq(0, end.time, by=1), out[,4], col="pink") # V2
# lines(seq(0, end.time, by=1), out[,5], col="violet") # V3
# 
# # dog quarantine
# plot(seq(0, end.time, by=1), rowSums(out[,8:11]), type="l") # all in qurantine
# lines(seq(0, end.time, by=1), rowSums(counts[,10:11]), col="red") # exp+inf in quarantine that will be removed
# 
# # dog infection in quarantine
# # plot(seq(0, end.time, by=1), counts$Qs, type="l", ylim=c(0,max(counts$Qs))) # Qs
# # lines(seq(0, end.time, by=1), counts$Qer, col="green") # Qer
# # lines(seq(0, end.time, by=1), counts$Qeb, col="blue") # Qeb
# # lines(seq(0, end.time, by=1), counts$Qi, col="darkblue") # Qi
# 
# plot(seq(0, end.time, by=1), out[,8], type="l", ylim=c(0,max(out[,8]))) # Qs
# lines(seq(0, end.time, by=1), out[,9], col="green") # Qer
# lines(seq(0, end.time, by=1), out[,10], col="blue") # Qeb
# lines(seq(0, end.time, by=1), out[,11], col="darkblue") # Qi
# 
# # dog infection in dogs and humans
# plot(seq(0, end.time, by=1), counts$I, type="l", ylim=c(0,max(counts$I))) # I dogs
# lines(seq(0, end.time, by=1), counts$Ih, col="red") # I humans
# 
# plot(seq(0, end.time, by=1), out[,7], type="l", ylim=c(0,max(out[,7]))) # I dogs
# lines(seq(0, end.time, by=1), out[,15], col="red") # I humans
# 
# plot(seq(0, end.time, by=1), out[,6], type="l", ylim=c(0,max(out[,6]))) # E dogs
# lines(seq(0, end.time, by=1), out[,7], col="red") # I dogs
# 
# plot(seq(0, end.time, by=1), counts$Rill, type="l", ylim=c(0,max(counts$Rill))) # dead dogs
# lines(seq(0, end.time, by=1), counts$Rh, col="red") # dead humans
# 
# plot(seq(0, end.time, by=1), counts$Rill, type="l") # Rabid dead dogs
# lines(seq(0, end.time, by=1), counts$I, col="red") # Rabid alive dogs
# sum(counts$I)
# sum(counts$Rill)
# 
# # plot(seq(0, end.time, by=1), out[,12], type="l") # dead dogs
# # lines(seq(0, end.time, by=1), out[,16], col="red") # dead humans
# 
# # PEP uptake
# plot(seq(0, end.time, by=1), counts$Vhs, type="l") # Vhs
# lines(seq(0, end.time, by=1), counts$Vhe, col="blue") # Vhe
# lines(seq(0, end.time, by=1), counts$Eh, col="red") # Eh (not given PEP)
# sum(counts$Eh)
# sum(counts$Ih)
# sum(counts$Rh)
# 
# 
# 
# 
