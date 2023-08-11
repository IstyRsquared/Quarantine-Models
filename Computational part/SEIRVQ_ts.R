## Quarantine model 
## Name: Isty Rysava
## Date: 11/08/2023
## Code: Reads in simulation outputs draws monthly ts of human and dog cases

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

## Libraries 
library(harrypotter)
library(colorspace)
library(tidyverse)

## Data
final_list_ts <- readRDS("output/MS_monthly_infection_ts.Rdata") 

### TIME SERIES (of dead dogs and humans)

# ## colours
# pal <- hp(n = 10, house = "Slytherin")
# colgreen <- pal[10]; colgreenlight <- pal[3]
# pal <- hp(n = 10, option = "LunaLovegood")
# colpink <- pal[7]; colpinkdark <- pal[10]
# 
# ## Plot only for R0=1.2 or 1.3
# Rsubset <- filter(params_grid, R0==1.3)
# for(idx in 1:nrow(Rsubset)){
#   R0 <- Rsubset[idx,][1]
#   sqc <- Rsubset[idx,][2]
#   vac <- Rsubset[idx,][3]
#   path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".pdf")
#   
#   ind <- which(paste(R0, sqc, vac)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
#   stats_mthly <- final_list_ts[[ind]]
#   dogs <- stats_mthly[[1]]
#   humans <- stats_mthly[[4]]
#   
#   pdf(path, width=6, height=4)
#   par(mar=c(5,5.5,1,1.5))
#   # main dogs
#   x<-1:nrow(dogs)
#   plot(dogs$mean~x,type="l", ylab="",xlab="", axes=F,ylim=c(0, 120), xaxt="n", 
#        bty="n", main="")
#   title(xlab = "Time (months)", line = 2, cex.lab=.7)            
#   title(ylab = "Canine rabies cases", line = 2, cex.lab=.7)   
#   axis(1, at=seq(1,61,12), labels=paste("Jan",c(2018, 2019, 2020, 2021, 2022, 2023)), cex.axis=0.7)
#   axis(2, cex.axis=0.7)
#   axis(1,at=seq(1,61,3),lab=rep("",length(seq(1,61,3))),tck=-0.01)
#   polygon(c(1:nrow(dogs),nrow(dogs):1),c(dogs$upperPI,rev(dogs$lowerPI)),col=colpink,border=F)
#   lines(dogs$mean~x,col=colpinkdark,lwd=2)
#   
#   # inset humans
#   par(new=TRUE, mar=c(0,0,0,0),mfrow=c(1,1), plt=c(0.67, 0.93, 0.67, 0.89), cex=0.72)
#   x<-1:nrow(humans)
#   plot(humans$mean~x,type="l",ylab="",xlab="",axes=F,ylim=c(0,max(humans$upperPI)), xaxt="n", bty="n")
#   title(ylab = "Human rabies cases", line = 2, cex.lab=.7)   
#   axis(1,at=seq(1,61,12),labels=c(2018, 2019, 2020, 2021, 2022, 2023), cex.axis=0.7, tck=-0.05)
#   axis(1, at=seq(1,61,3), lab=rep("",length(seq(1,61,3))),tck=-0.01)
#   axis(2, cex.axis=.7, tck=-0.05)
#   polygon(c(1:nrow(humans),nrow(humans):1),c(humans$upperPI,rev(humans$lowerPI)),col=colgreenlight,border=F)
#   lines(humans$mean~x,col=colgreen,lwd=2)
#   dev.off()
# }

## Plot all individually
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

path <- "figs/ts/MS_sim_runs_R01.3_sqcall_vcall.pdf"
pdf(path, width=5, height=8)
par(mfrow=c(3,1))

for(idx in 1:length(sqcs)){
  sqc <- sqcs[idx]
  # path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqc", sqc, "_vcall.pdf")
  # pdf(path, width=6, height=4)
  
  for(i in 1:length(vacs)){
    ind <- which(paste(R0, sqc, vacs[i])==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
    
    stats_mthly <- final_list_ts[[ind]]
    dogs <- stats_mthly[[1]]
    dogs <- dogs[(7:60),]
    # humans <- stats_mthly[[4]]
    
    # par(mar=c(5,5.5,1,1.5))
    if(i==1){
      x<-1:nrow(dogs)
      plot(dogs$mean~x,type="l", ylab="", xlab="", axes=F, ylim=c(0, 50), xaxt="n", bty="n")
      title(xlab = "Time", line = 2, cex.lab=.9)            
      title(ylab = "Monthly no. of infected dogs (E+I)", line = 2, cex.lab=.9)   
      axis(1, at=seq(7,55,12), labels=paste("Year", 2:6), cex.axis=0.8)
      axis(2, cex.axis=0.8)
      axis(1, at=seq(1,54,3), lab=rep("",length(seq(1,54,3))), tck=-0.01)
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI, rev(dogs$lowerPI)), col=borderq4[i], border=F)
      lines(dogs$mean~x, col=linesq4[i],lwd=1)
    }else{
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI, rev(dogs$lowerPI)), col=borderq4[i], border=F)
      lines(dogs$mean~x, col=linesq4[i],lwd=1)
    }
    
    # add legend
    if(idx==3 & i==4){
      legend(48, 50, title="Vaccination", legend=c("0%", "25%", "50%", "75%"),
             col=linesq4, lty=c(1,1,1,1), cex=0.7, box.lty=0)
    }
  }
  # dev.off()
}
dev.off()

path <- "figs/ts/MS_sim_runs_R01.3_sqcall_vcallB.pdf"
mains <- c("A/", "B/", "C/", "D/")
borderq4 <- c(alpha(q4[1:2], 0.1), alpha(q4[3], 0.3), alpha(q4[4], 0.3))
pdf(path, width=8, height=6)
par(mfrow=c(2,2))

for(idx in 1:length(vacs)){
  vacc <- vacs[idx]
  
  for(i in 1:length(sqcs)){
    ind <- which(paste(R0, sqcs[i], vacc)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
    
    stats_mthly <- final_list_ts[[ind]]
    dogs <- stats_mthly[[1]]
    dogs <- dogs[(7:60),]
    
    if(i==1){
      x<-1:nrow(dogs)
      plot(dogs$mean~x,type="l", ylab="", xlab="", axes=F, ylim=c(0, 50), xaxt="n", bty="n")
      title(xlab = "Time", line = 2, cex.lab=.8)     
      title(mains[idx], adj = 0, line = 1, cex.main=.8)
      title(ylab = "Monthly no. of infected dogs (E+I)", line = 2, cex.lab=.9)   
      axis(1, at=seq(7,55,12), labels=paste("Year", 2:6), cex.axis=0.8)
      axis(2, cex.axis=0.8)
      axis(1, at=seq(1,54,3), lab=rep("",length(seq(1,54,3))), tck=-0.01)
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI, rev(dogs$lowerPI)), col=borderq4[i], border=F)
      lines(dogs$mean~x, col=linesq4[i],lwd=1)
    }else{
      polygon(c(1:nrow(dogs), nrow(dogs):1), c(dogs$upperPI, rev(dogs$lowerPI)), col=borderq4[i], border=F)
      lines(dogs$mean~x, col=linesq4[i],lwd=1)
    }
    
    # add legend
    if(idx==4 & i==3){
      legend(35, 50, title="Quarantine scenario", legend=c("1", "2", "3"),
             col=linesq4, lty=c(1,1,1), cex=0.65, box.lty=0)
    }
  }
}
dev.off()

