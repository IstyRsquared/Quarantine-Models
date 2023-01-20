## Quarantine model 
## Name: Isty Rysava
## Date: 18/11//22
## Code: Reads in simulation outputs draws A) monthly ts of human and dog cases, 
## B) boxplot of monthly infection

rm(list=ls())

# setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

library(harrypotter)
library(tidyverse)

## Params set up
end.time <- 52*5
nsim <- 1000

R0s <- seq(1, 2, 0.1)
sqcs <- 1:3
vc.temp <-  c(0, 0.25, 0.5, 0.75)
params_grid <- expand.grid(list(R0 = R0s, # reproductive number
                                sqc = sqcs, # quarantine scenarios
                                vc = vc.temp)) # vaccination scenarios

final_frame_box <- c()
final_list_ts <- vector("list", length=nrow(params_grid))

weeks <- seq(as.Date("2018-01-01"), as.Date("2022-12-31"), by="week")
months <- ((as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$year-118)*12) + as.POSIXlt(strptime(weeks, format="%Y-%m-%d"))$mon+1

for(idx in 1:nrow(params_grid)){
  R0 <- params_grid[idx,][1]
  sqc <- params_grid[idx,][2]
  vac <- params_grid[idx,][3]
  my.files <- readRDS(paste0("output/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".Rdata"))
  
  ## Aggregate by month 
  sum_mthly <- vector("list", length(my.files))
  names(my.files) <- names(sum_mthly) <- list("deadD", "expD", "infD", "deadH")

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
  
  final_list_ts[[idx]] <- stats_mthly
  print(idx)
}

# saveRDS(final_list_ts, "output/MS_monthly_infection_ts.Rdata") 
# write.csv(final_frame_box, "output/MS_monthly_infection_boxplot.csv", row.names=F) 

#################################################################################################################################################
### TIME SERIES (of dead dogs and humans)
final_list_ts <- readRDS("output/MS_monthly_infection_ts.Rdata") 

## colours
pal <- hp(n = 10, house = "Slytherin")
colgreen <- pal[10]; colgreenlight <- pal[3]
pal <- hp(n = 10, option = "LunaLovegood")
colpink <- pal[7]; colpinkdark <- pal[10]

## Plot only for R0=1.2 or 1.3
Rsubset <- filter(params_grid, R0==1.3)
for(idx in 1:nrow(Rsubset)){
  R0 <- Rsubset[idx,][1]
  sqc <- Rsubset[idx,][2]
  vac <- Rsubset[idx,][3]
  path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".pdf")
  
  ind <- which(paste(R0, sqc, vac)==paste(params_grid$R0, params_grid$sqc, params_grid$vc))
  stats_mthly <- final_list_ts[[ind]]
  dogs <- stats_mthly[[1]]
  humans <- stats_mthly[[4]]
  
  pdf(path, width=6, height=4)
  par(mar=c(5,5.5,1,1.5))
  # main dogs
  x<-1:nrow(dogs)
  plot(dogs$mean~x,type="l", ylab="",xlab="", axes=F,ylim=c(0, max(dogs$upperPI)), xaxt="n", 
       bty="n", main=paste0("Vaccination ", vac*100, "%"), cex.main=.7)
  title(xlab = "Time (months)", line = 2, cex.lab=.7)            
  title(ylab = "Canine rabies cases", line = 2, cex.lab=.7)   
  # plot(dogs$mean~x,type="l",cex.lab=1,ylab="Canine rabies cases",xlab="Time (months)", axes=F,ylim=c(0,130), xaxt="n", bty="n")
  axis(1, at=seq(1,61,12), labels=paste("Jan",c(2018, 2019, 2020, 2021, 2022, 2023)), cex.axis=0.7)
  axis(2, cex.axis=0.7)
  axis(1,at=seq(1,61,3),lab=rep("",length(seq(1,61,3))),tck=-0.01)
  polygon(c(1:nrow(dogs),nrow(dogs):1),c(dogs$upperPI,rev(dogs$lowerPI)),col=colpink,border=F)
  lines(dogs$mean~x,col=colpinkdark,lwd=2)
  
  # inset humans
  par(new=TRUE, mar=c(0,0,0,0),mfrow=c(1,1), plt=c(0.67, 0.93, 0.67, 0.89), cex=0.72)
  x<-1:nrow(humans)
  plot(humans$mean~x,type="l",ylab="",xlab="",axes=F,ylim=c(0,max(humans$upperPI)), xaxt="n", bty="n")
  # plot(humans$mean~x,type="l",cex.lab=1,ylab="Human rabies cases",xlab="Time (months)",axes=F,ylim=c(0,10), xaxt="n", bty="n")
  # title(xlab = "Time (months)", line = 2)            
  title(ylab = "Human rabies cases", line = 2, cex.lab=.7)   
  axis(1,at=seq(1,61,12),labels=c(2018, 2019, 2020, 2021, 2022, 2023), cex.axis=0.7, tck=-0.05)
  axis(1, at=seq(1,61,3), lab=rep("",length(seq(1,61,3))),tck=-0.01)
  axis(2, cex.axis=.7, tck=-0.05)
  polygon(c(1:nrow(humans),nrow(humans):1),c(humans$upperPI,rev(humans$lowerPI)),col=colgreenlight,border=F)
  lines(humans$mean~x,col=colgreen,lwd=2)
  dev.off()
}

## Plot all individually
for(idx in 1:nrow(params_grid)){
  R0 <- params_grid[idx,][1]
  sqc <- params_grid[idx,][2]
  vac <- params_grid[idx,][3]
  path <- paste0("figs/ts/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".pdf")
  
  stats_mthly <- final_list_ts[[idx]]
  dogs <- stats_mthly[[1]]
  humans <- stats_mthly[[4]]
  
  pdf(path, width=6, height=4)
  par(mar=c(5,5.5,1,1.5), cex=0.8)
  # main dogs
  x<-1:nrow(dogs)
  plot(dogs$mean~x,type="l",cex.lab=1,ylab="",xlab="",axes=F,ylim=c(0,max(dogs$upperPI)), xaxt="n", bty="n",
       main=paste0("R0=", params_grid[idx,1], " ", params_grid[idx,2], " ", "vc=", params_grid[idx,3]), cex.main=.8)
  title(xlab = "Time (months)", line = 2)            
  title(ylab = "Canine rabies cases", line = 2)   
  # plot(dogs$mean~x,type="l",cex.lab=1,ylab="Canine rabies cases",xlab="Time (months)", axes=F,ylim=c(0,130), xaxt="n", bty="n")
  axis(1, at=seq(1,61,12), labels=paste("Jan",c(2018, 2019, 2020, 2021, 2022, 2023)), cex.axis=0.9)
  axis(2, cex.axis=0.9)
  axis(1,at=seq(1,61,3),lab=rep("",length(seq(1,61,3))),tck=-0.01)
  polygon(c(1:nrow(dogs),nrow(dogs):1),c(dogs$upperPI,rev(dogs$lowerPI)),col=colpink,border=F)
  lines(dogs$mean~x,col=colpinkdark,lwd=2)
  
  # inset humans
  par(new=TRUE, mar=c(0,0,0,0),mfrow=c(1,1), plt=c(0.67, 0.93, 0.67, 0.89), cex=0.72)
  x<-1:nrow(humans)
  plot(humans$mean~x,type="l",cex.lab=1,ylab="",xlab="",axes=F,ylim=c(0,max(humans$upperPI)), xaxt="n", bty="n")
  #plot(humans$mean~x,type="l",cex.lab=1,ylab="Human rabies cases",xlab="Time (months)",axes=F,ylim=c(0,10), xaxt="n", bty="n")
  title(xlab = "Time (months)", line = 2)            
  title(ylab = "Human rabies cases", line = 2)   
  axis(1,at=seq(1,61,12),labels=c(2018, 2019, 2020, 2021, 2022, 2023), cex.axis=0.8)
  axis(1, at=seq(1,61,3), lab=rep("",length(seq(1,61,3))),tck=-0.01)
  axis(2, cex.axis=.8)
  polygon(c(1:nrow(humans),nrow(humans):1),c(humans$upperPI,rev(humans$lowerPI)),col=colgreenlight,border=F)
  lines(humans$mean~x,col=colgreen,lwd=2)
  dev.off()
}

## Calculate human cases
# low: 4, high: 6
idx=6
stats_mthly <- final_list_ts[[idx]]
dogs <- stats_mthly[[1]]
humans <- stats_mthly[[4]]
sum(humans$mean) 
# R0=1.3, low: 65.076, high: 47.129
# (65.076 - 47.129) / (65.076/100) # 28%
# R0=1.2, low: 60.555, high: 45.101
# (60.555 - 45.101) / (60.555/100) # 26%

sum(dogs$mean) 
# R0=1.3, low: 958.651, high: 817.13
# (958.651 - 817.13) / (958.651/100) # 15%
# R0=1.2, low: 886.805, high: 776.412
# (886.805 - 776.412) / (886.805/100) # 13%

### BOXPLOT (of exposed and infectious dogs)
final_frame_box <- read.csv("output/MS_monthly_infection_boxplot.csv")

## colours
palbox <- hp(n = 4, house = "RonWeasley")
plot(1:4, 1:4, col=palbox, pch=16, cex=3)

library("colorspace")
q4 <- sequential_hcl(4, "PuBuGn")
plot(1:4, 1:4, col=q4, pch=16, cex=3)

mygray <- alpha("gray", 0.1)

## Set up theme
theme_set(theme_bw() +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.border=element_blank(),
                  axis.text=element_text(size=5),
                  axis.title.x = element_text(size = 5),
                  axis.title.y = element_text(size = 5),
                  # legend.title=element_blank(),
                  legend.text = element_text(size = 4),
                  text=element_text(size=5),
                  legend.key.size = unit(0.3, 'cm')))

# data <- data.frame(team=rep(c('A', 'B', 'C'), each=50),
#                    program=rep(c('low', 'high'), each=25),
#                    values=seq(1:150)+sample(1:100, 150, replace=TRUE))
# tibble(data)
# ggplot(data, aes(x=team, y=values, fill=program)) +
#   geom_boxplot()

# labs
supp.labs <- paste("R0 =", unique(final_frame_box$R0))
to_string <- as_labeller(c(`1` = supp.labs[1], `1.1` = supp.labs[2], `1.2` = supp.labs[3], `1.3` = supp.labs[4], `1.4` = supp.labs[5],
                           `1.5` = supp.labs[6], `1.6` = supp.labs[7], `1.7` = supp.labs[8], `1.8` = supp.labs[9], `1.9` = supp.labs[10], 
                           `2` = supp.labs[11]))
## Draw plot
tibble(final_frame_box)
final_frame_box <- tibble(final_frame_box) %>%
  mutate(R0 = as.factor(R0)) %>%
  mutate(Vaccination = as.factor(Vaccination)) %>%
  mutate(Quarantine = as.factor(Quarantine))

# test <- filter(final_frame_box, R0==1.3)
# ggplot(test, aes(x=Quarantine, y=expD+infD, fill=Vaccination)) +
#   geom_boxplot()
p_box <- ggplot(data = final_frame_box, aes(x=Quarantine, y=expD+infD, fill=Vaccination)) +
  geom_boxplot(outlier.size = 0.05, lwd=0.1, outlier.color=mygray) +
  # geom_boxplot(outlier.shape = NA) +
  facet_wrap(~R0, labeller = to_string, scales="free") +
  scale_fill_manual(name="Vaccination", values=c(q4), labels=c("0%", "25%", "50%", "75%")) + 
  ylab("Monthly no. of infected dogs (E+I)") + xlab("Quarantine scenario") +
  theme(strip.background = element_blank())

g1 <- p_box %+% dplyr::filter(final_frame_box, R0 == 1.3) + theme(legend.position = "none")
g2 <- p_box %+% dplyr::filter(final_frame_box, R0 != 1.3) + facet_wrap(~R0, nrow=2, scales="free")

Infection_CompGrid <- gridExtra::grid.arrange(g1, g2,
                                          layout_matrix = 
                                            matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                                     1, 1, 1, 2, 2, 2, 2, 2, 2),
                                                   byrow = TRUE, nrow = 2))

ggsave("figs/Infection_CompGrid.png", plot=Infection_CompGrid, width = 20, height = 9, units = "cm")









