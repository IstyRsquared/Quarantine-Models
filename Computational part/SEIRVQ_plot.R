## Hong KONG model development
## Name: Isty Rysava
## Date: 27/06/18
## Code: reads in simulation outputs of 3 quarantine scenarios and draws monthly ts of human and dog cases

rm(list=ls())
setwd("~/Dropbox/HongKong")
library(harrypotter)

### Data 
sc2 <- readRDS("counts_thesis_scenario2extravacc.Rdata")

### Set up a desired simulation to process ##
sc <- sc2
name <- "sc2counts"
nsim <- length(sc)
time <- 1:length(sc[[1]][,1])
treatment <- "extravacc" # original, halfvacc, extravacc

pathA <- paste0("output/", name, "_")
pathB <- paste0("output/", name, "_monthly_")
pathC <- paste("output/", name, "_statsmonthly_")
pathD <- paste0("figs/monthly_ts_", name, "_", treatment, ".pdf")

## Substract populations of interest ##
# initialize matrices
deadD <- matrix(NA, ncol=nsim, nrow=length(time))
deadH <- matrix(NA, ncol=nsim, nrow=length(time))

# populate
for(i in 1:nsim){
  colnames(sc[[i]]) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh",
                         "Ih", "Rh", "Vhs", "Vhe")
  deadD[,i] <- sc[[i]][,"Rill"]
  deadH[,i] <- sc[[i]][,"Rh"]
}

# save
# write.csv(deadD, paste0(pathA, "deadD.csv"), row.names=F)
# write.csv(deadH, paste0(pathA, "deadH.csv"), row.names=F)

## Aggregate by month ##
weeks <- seq(as.Date("2018-01-01"), as.Date("2022-12-31"), by="week")
months <- ((as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$year-118)*12) + as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$mon+1

# initialize lists
my.files <- list(deadD, deadH)
sum_mthly <- vector("list", length(my.files))
names(my.files) <- list("deadD", "deadH")
names(sum_mthly) <- list("deadD", "deadH")

# aggreagte
for(i in 1: length(my.files)){
  df.temp <- data.frame(my.files[[i]])
  sum_mthly[[i]] <- rowsum(df.temp, group=months)
}

# save 
# write.csv(sum_mthly[[4]], paste0(pathB, "deadD.csv"), row.names=F)
# write.csv(sum_mthly[[5]], paste0(pathB, "deadH.csv"), row.names=F)

### Calculate statistics ###
# initialize list
stats_mthly <- vector("list", length(sum_mthly))
names(stats_mthly) <- list("deadD", "deadH")

for(i in 1:length(sum_mthly)){
  df.temp <- data.frame(sum_mthly[[i]])
  monthly_stats <- data.frame(mean=as.numeric(rowMeans(df.temp)))
  for(j in 1:nrow(monthly_stats)){
    monthly_stats$upperPI[j] <- as.numeric(as.character(sort(df.temp[j,])[round(0.975*nsim)]))
    monthly_stats$lowerPI[j] <- as.numeric(as.character(sort(df.temp[j,])[round(0.275*nsim)]))
  }
  stats_mthly[[i]] <- monthly_stats
}

# save 
# write.csv(stats_mthly[[1]], paste0(pathC, "deadD.csv"), row.names=F)
# write.csv(stats_mthly[[2]], paste0(pathC, "deadH.csv"), row.names=F)

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










