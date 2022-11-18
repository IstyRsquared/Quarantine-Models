## Quarantine model 
## Name: Isty Rysava
## Date: 18/11//22
## Code: Reads in simulation outputs draws A) monthly ts of human and dog cases, 
## B) boxplot of monthly infection

rm(list=ls())

setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
library(harrypotter)
library(tidyverse)

### Data 
allout <- readRDS("output/MS_sim_runs_test.Rdata")
final_frame_box <- c()
final_frame_ts <- c()

idx <- c(37, 38)

## Params set up
end.time <- 52*5
nsim <- 1000

R0s <- seq(1, 2, 0.1)
sqcs <- paste("Scenario", 1:3)
vc.temp <-  c("0%", "25%", "50%", "75%")
params_grid <- expand.grid(list(R0 = R0s, # reproductive number
                                sqc = sqcs, # quarantine scenarios
                                vc = vc.temp)) # vaccination scenarios

for(idx in 1:nrow(params_grid)){
  out.temp <- allout[[idx]]
  
  ## Aggregate by month ##
  weeks <- seq(as.Date("2018-01-01"), as.Date("2022-12-31"), by="week")
  months <- ((as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$year-118)*12) + as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$mon+1
  
  my.files <- out.temp
  sum_mthly <- vector("list", length(my.files))
  names(my.files) <- list("deadD", "expD", "infD", "deadH")
  names(sum_mthly) <- list("deadD", "expD", "infD", "deadH")
  
  for(i in 1: length(my.files)){
    df.temp <- data.frame(my.files[[i]])
    sum_mthly[[i]] <- rowsum(df.temp, group=months)
  }
  
  ### Boxplot figs (keep as is)
  final_frame_box <- rbind(final_frame_box, 
                       data.frame(R0=params_grid$R0[idx], Quarantine=params_grid$sqc[idx], Vaccination=params_grid$vc[idx],
                                  deadD=as.vector(t(sum_mthly[[1]])), expD=as.vector(t(sum_mthly[[2]])), 
                                  infD=as.vector(t(sum_mthly[[3]])), deadH=as.vector(t(sum_mthly[[4]]))))
  
  ### Time series figs (average over sims)
  stats_mthly <- vector("list", length(sum_mthly))
  names(stats_mthly) <- list("deadD", "expD", "infD", "deadH")
  
  for(i in 1:length(sum_mthly)){
    df.temp <- data.frame(sum_mthly[[i]])
    monthly_stats <- data.frame(mean=as.numeric(rowMeans(df.temp)))
    for(j in 1:nrow(monthly_stats)){
      monthly_stats$upperPI[j] <- as.numeric(as.character(sort(df.temp[j,])[round(0.975*nsim)]))
      monthly_stats$lowerPI[j] <- as.numeric(as.character(sort(df.temp[j,])[round(0.275*nsim)]))
    }
    stats_mthly[[i]] <- monthly_stats
  }
  
  ## HERE!!!
  stats_mthly[[1]]
  final_frame_ts <- rbind(final_frame_ts, 
                           data.frame(R0=params_grid$R0[idx], Quarantine=params_grid$sqc[idx], Vaccination=params_grid$vc[idx],
                                      deadD=as.vector(t(sum_mthly[[1]])), expD=as.vector(t(sum_mthly[[2]])), 
                                      infD=as.vector(t(sum_mthly[[3]])), deadH=as.vector(t(sum_mthly[[4]]))))
  
  
  
}

#################################################################################################################################################
data <- data.frame(team=rep(c('A', 'B', 'C'), each=50),
                   program=rep(c('low', 'high'), each=25),
                   values=seq(1:150)+sample(1:100, 150, replace=TRUE))
head(data)
ggplot(data, aes(x=team, y=values, fill=program)) + 
  geom_boxplot() 



### Prepare colours
pal <- hp(n = 10, house = "Slytherin")
plot(1:10, 1:10, col=pal, pch=16, cex=3)
colgreen <- pal[10]
colgreenlight <- pal[3]

pal <- hp(n = 10, option = "LunaLovegood")
plot(1:10, 1:10, col=pal, pch=16, cex=3)
colpink <- pal[7]
colpinkdark <- pal[10]

### Plot time series ###
# of dead dogs and humans
dogs <- stats_mthly[[1]]
colSums(dogs)
humans <- stats_mthly[[2]]
colSums(humans)

pdf(pathD, width=6, height=4)
par(mar=c(5,5.5,1,1.5), cex=0.9)
# main dogs
x<-1:nrow(dogs)
#plot(dogs$mean~x,type="l",cex.lab=1.2,ylab="Canine rabies cases",xlab="Time (months)",axes=F,ylim=c(0,max(dogs$upperPI)), xaxt="n", bty="n")
plot(dogs$mean~x,type="l",cex.lab=1,ylab="Canine rabies cases",xlab="Time (months)", axes=F,ylim=c(0,130), xaxt="n", bty="n")
axis(1, at=seq(1,61,12), labels=paste("Jan",c(2018, 2019, 2020, 2021, 2022, 2023)), cex.axis=0.9)
axis(2, cex.axis=0.9)
axis(1,at=seq(1,61,3),lab=rep("",length(seq(1,61,3))),tck=-0.01)
polygon(c(1:nrow(dogs),nrow(dogs):1),c(dogs$upperPI,rev(dogs$lowerPI)),col=colpink,border=F)
lines(dogs$mean~x,col=colpinkdark,lwd=2)

# inset humans
par(new=TRUE, mar=c(0,0,0,0),mfrow=c(1,1), plt=c(0.67, 0.93, 0.67, 0.89), cex=0.72)
x<-1:nrow(humans)
#plot(humans$mean~x,type="l",cex.lab=1.2,ylab="Human rabies cases",xlab="Time (months)",axes=F,ylim=c(0,max(humans$upperPI)), xaxt="n", bty="n")
plot(humans$mean~x,type="l",cex.lab=1,ylab="Human rabies cases",xlab="Time (months)",axes=F,ylim=c(0,10), xaxt="n", bty="n")
axis(1,at=seq(1,61,12),labels=c(2018, 2019, 2020, 2021, 2022, 2023), cex.axis=0.8)
axis(1, at=seq(1,61,3), lab=rep("",length(seq(1,61,3))),tck=-0.01)
axis(2, cex.axis=.8)
polygon(c(1:nrow(humans),nrow(humans):1),c(humans$upperPI,rev(humans$lowerPI)),col=colgreenlight,border=F)
lines(humans$mean~x,col=colgreen,lwd=2)
dev.off()










