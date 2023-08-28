## Quarantine model 
## Name: Isty Rysava
## Date: 11/08/2023
## Code: Reads in simulation outputs and draws monthly ts of human and dog cases

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

## Libraries 
library(harrypotter)
library(colorspace)
library(tidyverse)

## Params set up
end.time <- 52*10
nsim <- 1000

R0s <- seq(1, 2, 0.1)
sqcs <- 1:3 # 1=0, 2=1.5, 3=3
vc.temp <-  c(0, 0.25, 0.5, 0.75)
params_grid <- expand.grid(list(R0 = R0s, # reproductive number
                                sqc = sqcs, # incursion scenarios
                                vc = vc.temp)) # vaccination scenarios

final_frame_box <- c()
final_frame_box_burnout <- c() # remove burn-in period (1st 6 months)
final_list_ts <- vector("list", length=nrow(params_grid))

weeks <- seq(as.Date("2018-01-01"), as.Date("2027-12-31"), by="week")[2:522]
months <- ((as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$year-118)*12) + as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$mon+1

for(idx in 1:nrow(params_grid)){
  R0 <- params_grid[idx,][1]
  sqc <- params_grid[idx,][2]
  vac <- params_grid[idx,][3]
  my.files <- readRDS(paste0("output/incs/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".Rdata"))
  
  ## Aggregate by month 
  sum_mthly <- vector("list", length(my.files))
  names(my.files) <- names(sum_mthly) <- list("deadD", "expD", "infD", "deadH", "pop")
  
  for(i in 1: length(my.files)){
    df.temp <- data.frame(my.files[[i]])
    sum_mthly[[i]] <- rowsum(df.temp, group=months)
  }
  
  ### Boxplot figs (keep as is)
  # with burn-in
  final_frame_box <- rbind(final_frame_box, 
                           data.frame(R0=params_grid$R0[idx], Quarantine=params_grid$sqc[idx], Vaccination=params_grid$vc[idx],
                                      deadD=as.vector(t(sum_mthly[[1]])), expD=as.vector(t(sum_mthly[[2]])), 
                                      infD=as.vector(t(sum_mthly[[3]])), deadH=as.vector(t(sum_mthly[[4]])), 
                                      pop=as.vector(t(sum_mthly[[5]]))))
  
  
  # without burn-in
  deadD_df <- t(sum_mthly[[1]])
  expD_df <- t(sum_mthly[[2]])
  infD_df <- t(sum_mthly[[3]])
  deadH_df <- t(sum_mthly[[4]])
  pop_df <- t(sum_mthly[[5]])
  final_frame_box_burnout <- rbind(final_frame_box_burnout, 
                                   data.frame(R0=params_grid$R0[idx], Quarantine=params_grid$sqc[idx], Vaccination=params_grid$vc[idx],
                                              deadD=as.vector(deadD_df[,7:60]), expD=as.vector(expD_df[,7:60]), 
                                              infD=as.vector(infD_df[,7:60]), deadH=as.vector(deadH_df[,7:60]), 
                                              pop=as.vector(pop_df[,7:60])))
  
  ### Time series figs (average over sims)
  stats_mthly <- vector("list", length(sum_mthly))
  names(stats_mthly) <- list("deadD", "expD", "infD", "deadH", "pop")
  
  for(i in 1:length(sum_mthly)){
    df.temp <- data.frame(sum_mthly[[i]])
    monthly_stats <- data.frame(mean=as.numeric(rowMeans(df.temp)))
    for(j in 1:nrow(monthly_stats)){
      monthly_stats$upperPI[j] <- as.numeric(as.character(sort(as.numeric(df.temp[j,]))[round(0.975*nsim)]))
      monthly_stats$lowerPI[j] <- as.numeric(as.character(sort(as.numeric(df.temp[j,]))[round(0.275*nsim)]))
    }
    stats_mthly[[i]] <- monthly_stats
  }
  
  final_list_ts[[idx]] <- stats_mthly
  print(idx)
}

# saveRDS(final_list_ts, "output/incs/MS_monthly_infection_tsinc_10yrs.Rdata")
# colnames(final_frame_box)[2] <- "Incursion"
# write.csv(final_frame_box_burnout, "output/incs/MS_monthly_infection_boxplotinc_burnout_10yrs.csv", row.names=F) 

#################################################################################################################################################
### TIME SERIES 1: dead dogs across incursion rates
## colours
q4 <- sequential_hcl(4, "PuBuGn")
plot(1:4, 1:4, col=q4, pch=16, cex=3)
linesq4 <- c(q4[1:3], "pink")
plot(1:4, 1:4, col=linesq4, pch=16, cex=3)
borderq4 <- c(alpha(q4[1:2], 0.1), alpha(q4[3], 0.2), alpha(q4[4], 0.3))
plot(1:4, 1:4, col=borderq4, pch=16, cex=3)

R0 <- 1.3
sqcs <- 1:3
vacs <- c(0, 0.25, 0.5, 0.75)

final_list_ts <- readRDS("output/incs/MS_monthly_infection_tsinc_10yrs.Rdata")  

path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqcall_vcall_INCS.pdf")
mains <- c("A/", "B/", "C/", "D/")
borderq4 <- c(alpha(q4[1:2], 0.1), alpha(q4[3], 0.3), alpha(q4[4], 0.3))
pdf(path, width=8, height=6)
par(mfrow=c(2,2))

for(idx in 1:length(vacs)){
  vacc <- vacs[idx]
  
  for(i in 1:length(sqcs)){
    ind <- which(paste(R0, sqcs[i], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
    max.temp <- final_list_ts[[which(paste(R0, sqcs[3], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))]]
    max <- ceiling(max(max.temp[[1]]$upperPI/max.temp[[5]]$upperPI*10000, na.rm=T))
    
    stats_mthly <- final_list_ts[[ind]]
    dogs <- stats_mthly[[1]]
    pop <- stats_mthly[[5]]
    dogs <- dogs[(7:120),]
    pop <- pop[(7:120),]
    
    if(i==1){
      x<-1:nrow(dogs)
      plot(dogs$mean/pop$mean*10000~x,type="l", ylab="", xlab="", axes=F, ylim=c(0, max), xaxt="n", bty="n")
      title(xlab = "Time", line = 2, cex.lab=.8)     
      title(mains[idx], adj = 0, line = 1, cex.main=.8)
      title(ylab = "Monthly incidence of dead dogs / 10,000", line = 2, cex.lab=.9)   
      axis(1, at=seq(7,120,12), labels=paste("Year", 2:11), cex.axis=0.8) 
      axis(2, cex.axis=0.8)
      axis(1, at=seq(1,54,3), lab=rep("",length(seq(1,54,3))), tck=-0.01)
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI/pop$upperPI*10000, rev(dogs$lowerPI/pop$lowerPI*10000)), 
              col=borderq4[i], border=F)
      lines(dogs$mean/pop$mean*10000, col=linesq4[i],lwd=1)
    }else{
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI/pop$upperPI*10000, rev(dogs$lowerPI/pop$lowerPI*10000)), 
              col=borderq4[i], border=F)
      lines(dogs$mean/pop$mean*10000, col=linesq4[i],lwd=1)
    }
    
    # add legend
    if(idx==4 & i==3){
      legend(95, 1, title="Incursion rate", legend=c("1", "1.5", "3"),
             col=linesq4, lty=c(1,1,1), cex=0.65, box.lty=0)
    }
  }
}
dev.off()

#################################################################################################################################################
### TIME SERIES 2: incidence per 100,000 across incursion rates
final_list_ts <- readRDS("output/incs/MS_monthly_infection_tsinc_10yrs.Rdata")  

path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqcall_vcall_INCSincidence.pdf")
mains <- c("A/", "B/", "C/", "D/")
borderq4 <- c(alpha(q4[1:2], 0.1), alpha(q4[3], 0.3), alpha(q4[4], 0.3))
pdf(path, width=8, height=6)
par(mfrow=c(2,2))

for(idx in 1:length(vacs)){
  vacc <- vacs[idx]
  
  for(i in 1:length(sqcs)){
    ind <- which(paste(R0, sqcs[i], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
    max.temp <- final_list_ts[[which(paste(R0, sqcs[3], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))]]
    max <- ceiling(max((max.temp[[2]]$upperPI+max.temp[[3]]$upperPI)/max.temp[[5]]$upperPI*10000, na.rm=T))
    
    stats_mthly <- final_list_ts[[ind]]
    Edogs <- stats_mthly[[2]]; Edogs <- Edogs[(7:120),]
    Idogs <- stats_mthly[[3]]; Idogs <- Idogs[(7:120),]
    pop <- stats_mthly[[5]]; pop <- pop[(7:120),]
    dogs <- Edogs+Idogs
    
    if(i==1){
      x<-1:nrow(dogs)
      plot(dogs$mean/pop$mean*10000~x,type="l", ylab="", xlab="", axes=F, ylim=c(0, max), xaxt="n", bty="n")
      title(xlab = "Time", line = 2, cex.lab=.8)     
      title(mains[idx], adj = 0, line = 1, cex.main=.8)
      title(ylab = "Monthly incidence of infected dogs (E+I) / 10,000", line = 2, cex.lab=.9)   
      axis(1, at=seq(7,120,12), labels=paste("Year", 2:11), cex.axis=0.9) 
      axis(2, cex.axis=0.9)
      axis(1, at=seq(1,54,3), lab=rep("",length(seq(1,54,3))), tck=-0.01)
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI/pop$upperPI*10000, rev(dogs$lowerPI/pop$lowerPI*10000)), 
              col=borderq4[i], border=F)
      lines(dogs$mean/pop$mean*10000, col=linesq4[i],lwd=1)
    }else{
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI/pop$upperPI*10000, rev(dogs$lowerPI/pop$lowerPI*10000)), 
              col=borderq4[i], border=F)
      lines(dogs$mean/pop$mean*10000, col=linesq4[i],lwd=1)
    }
    
    # add legend
    if(idx==4 & i==3){
      legend(80, 3, title="Incursion scenario", legend=c("1", "2", "3"),
             col=linesq4, lty=c(1,1,1), lwd=c(2,2,2), cex=0.8, box.lty=0)
    }
  }
}
dev.off()

#################################################################################################################################################
### TIME SERIES 3: dead humans across incursion rates
final_list_ts <- readRDS("output/MS_monthly_infection_ts_10yrs.Rdata")  

path <- "figs/ts/MS_sim_runs_R01.3_sqcall_vcallHUMANS.pdf"
mains <- c("A/", "B/", "C/", "D/")
borderq4 <- c(alpha(q4[1:2], 0.1), alpha(q4[3], 0.3), alpha(q4[4], 0.3))
pdf(path, width=8, height=6)
par(mfrow=c(2,2))

for(idx in 1:length(vacs)){
  vacc <- vacs[idx]
  
  for(i in 1:length(sqcs)){
    ind <- which(paste(R0, sqcs[i], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
    
    stats_mthly <- final_list_ts[[ind]]
    humans <- stats_mthly[[4]]
    humans <- humans[(7:120),]

    if(i==1){
      x<-1:nrow(humans)
      plot(humans$mean~x,type="l", ylab="", xlab="", axes=F, ylim=c(0, 8), xaxt="n", bty="n")
      title(xlab = "Time", line = 2, cex.lab=.8)     
      title(mains[idx], adj = 0, line = 1, cex.main=.8)
      title(ylab = "Monthly human deaths due to rabies", line = 2, cex.lab=.9)   
      axis(1, at=seq(7,120,12), labels=paste("Year", 2:11), cex.axis=0.8) 
      axis(2, cex.axis=0.8)
      axis(1, at=seq(1,54,3), lab=rep("",length(seq(1,54,3))), tck=-0.01)
      polygon(c(1:nrow(humans), nrow(humans):1), c(humans$upperPI, rev(humans$lowerPI)), 
              col=borderq4[i], border=F)
      lines(humans$mean, col=linesq4[i],lwd=1)
    }else{
      polygon(c(1:nrow(humans), nrow(humans):1), c(humans$upperPI, rev(humans$lowerPI)), 
              col=borderq4[i], border=F)
      lines(humans$mean, col=linesq4[i],lwd=1)
    }
    
    # add legend
    if(idx==4 & i==3){
      legend(95, 8, title="Quarantine scenario", legend=c("1", "1.5", "3"),
             col=linesq4, lty=c(1,1,1), cex=0.65, box.lty=0)
    }
  }
}
dev.off()

