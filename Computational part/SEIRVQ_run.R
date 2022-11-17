## Quarantine model 
## Name: Isty Rysava
## Date: 08/02/18
## Code: runs  SEIR.tauleap model to test intervention scenarios (as per "thesis_params.R")

rm(list=ls())
setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

### Libraries 
source("R/TauLeap_SEIRfc_vacc_explicit_updated.R")
source("R/thesis_params.R")
source("R/thesis_vaccination_params.R")

### Parameters 
parameters <- vector("list", length=2)
# vaccweekly5yrs <- 0.5*vaccweekly5yrs # reduce vaccination coverage by 50% 
vaccweekly5yrs <- vaccweekly5yrs + (0.5*vaccweekly5yrs) # increase vaccination coverage by 50% 
parameters[[1]] <- c(0, vaccweekly5yrs) 
parameters[[2]] <- params3 # params1, params2, params3 (different quarantine treatments)

### Initial conditions 
initials <- c(S=S, V1=V1, V2=V2, V3=V3, E=E, I=I, Qs=0, Qer=0, Qeb=0, Qi=0,
              Rill=0, Sh=1324935, Eh=0, Ih=0,Rh=0, Vhs=0, Vhe=0) 
# humpop Google: 1,314,826 in 2015, mine: 1,324,935

end.time <- 52*5
nsim <- 1000
allout <- vector("list", length=nsim)
allcounts <- vector("list", length=nsim)

### Run the model: simulations
for(sim in 1:nsim){
  sierTL.out <- SEIR.tauleap(init = initials, pars = parameters,
                             end.time = end.time, tau=1) # 1/52 weekly steps
  allout[[sim]] <- as.matrix(sierTL.out$results)
  allcounts[[sim]] <- as.matrix(sierTL.out$counts)
  print(sim)
}

## original
# saveRDS(allout, "output/out_thesis_scenario1.Rdata") # params 1
# saveRDS(allcounts, "output/counts_thesis_scenario1.Rdata") # params 1
# saveRDS(allout, "output/out_thesis_scenario2.Rdata") # params 2
# saveRDS(allcounts, "output/counts_thesis_scenario2.Rdata") # params 2
# saveRDS(allout, "output/out_thesis_scenario3.Rdata") # params 3
# saveRDS(allcounts, "output/counts_thesis_scenario3.Rdata") # params 3

## 50% vaccination
# saveRDS(allout, "output/out_thesis_scenario1halfvacc.Rdata") # params 1
# saveRDS(allcounts, "output/counts_thesis_scenario1halfvacc.Rdata") # params 1
# saveRDS(allout, "output/out_thesis_scenario2halfvacc.Rdata") # params 2
# saveRDS(allcounts, "output/counts_thesis_scenario2halfvacc.Rdata") # params 2
# saveRDS(allout, "output/out_thesis_scenario3halfvacc.Rdata") # params 3
# saveRDS(allcounts, "output/counts_thesis_scenario3halfvacc.Rdata") # params 3

## extra 50% vaccination
# saveRDS(allout, "output/out_thesis_scenario1extravacc.Rdata") # params 1
# saveRDS(allcounts, "output/counts_thesis_scenario1extravacc.Rdata") # params 1
# saveRDS(allout, "output/out_thesis_scenario2extravacc.Rdata") # params 2
# saveRDS(allcounts, "output/counts_thesis_scenario2extravacc.Rdata") # params 2
# saveRDS(allout, "output/out_thesis_scenario3extravacc.Rdata") # params 3
# saveRDS(allcounts, "output/counts_thesis_scenario3extravacc.Rdata") # params 3

### Run the model: testing
end.time <- 52*5
sierTL.out <- SEIR.tauleap(init = initials, pars = parameters,
                           end.time = end.time, tau=1) # 1/52 weekly steps
out <- sierTL.out$results
counts <- sierTL.out$counts

out <- as.matrix(out)
counts <- as.data.frame(counts)
colnames(counts) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh",
                   "Ih", "Rh", "Vhs", "Vhe")
head(counts); tail(counts)

### Plot
# dog population
plot(seq(0, end.time, by=1), out[,2], type="l", ylim=c(0,max(out[,2]))) # sus
lines(seq(0, end.time, by=1), rowSums(out[,3:5]), col="orange") # all vacc

# dog vaccination
plot(seq(0, end.time, by=1), out[,3],  type="l", col="red", ylim=c(0,max(out[,3]))) # V1
lines(seq(0, end.time, by=1), out[,4], col="pink") # V2
lines(seq(0, end.time, by=1), out[,5], col="violet") # V3

# dog quarantine
plot(seq(0, end.time, by=1), rowSums(out[,8:11]), type="l") # all in qurantine
lines(seq(0, end.time, by=1), rowSums(counts[,10:11]), col="red") # exp+inf in quarantine that will be removed

# dog infection in quarantine
# plot(seq(0, end.time, by=1), counts$Qs, type="l", ylim=c(0,max(counts$Qs))) # Qs
# lines(seq(0, end.time, by=1), counts$Qer, col="green") # Qer
# lines(seq(0, end.time, by=1), counts$Qeb, col="blue") # Qeb
# lines(seq(0, end.time, by=1), counts$Qi, col="darkblue") # Qi

plot(seq(0, end.time, by=1), out[,8], type="l", ylim=c(0,max(out[,8]))) # Qs
lines(seq(0, end.time, by=1), out[,9], col="green") # Qer
lines(seq(0, end.time, by=1), out[,10], col="blue") # Qeb
lines(seq(0, end.time, by=1), out[,11], col="darkblue") # Qi

# dog infection in dogs and humans 
plot(seq(0, end.time, by=1), counts$I, type="l") # I dogs
lines(seq(0, end.time, by=1), counts$Ih, col="red") # I humans

plot(seq(0, end.time, by=1), out[,7], type="l") # I dogs
lines(seq(0, end.time, by=1), out[,15], col="red") # I humans

plot(seq(0, end.time, by=1), counts$Rill, type="l", ylim=c(0,max(counts$Rill))) # dead dogs
lines(seq(0, end.time, by=1), counts$Rh, col="red") # dead humans

plot(seq(0, end.time, by=1), counts$Rill, type="l") # Rabid dead dogs
lines(seq(0, end.time, by=1), counts$I, col="red") # Rabid alive dogs
sum(counts$I)
sum(counts$Rill)

# plot(seq(0, end.time, by=1), out[,12], type="l") # dead dogs
# lines(seq(0, end.time, by=1), out[,16], col="red") # dead humans: this is a large number, but probably because of so many infected dogs

# PEP uptake
plot(seq(0, end.time, by=1), counts$Vhs, type="l") # Vhs
lines(seq(0, end.time, by=1), counts$Vhe, col="blue") # Vhe
lines(seq(0, end.time, by=1), counts$Eh, col="red") # Eh (not given PEP)

sum(counts$Eh)
sum(counts$Ih)
sum(counts$Rh)




