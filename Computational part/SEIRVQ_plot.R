## Quarantine model 
## Name: Isty Rysava
## Date: 18/11/2022
## Code: Reads in simulation outputs draws boxplot of monthly infection 

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
sqcs <- 1:3
vc.temp <-  c(0, 0.25, 0.5, 0.75)
params_grid <- expand.grid(list(R0 = R0s, # reproductive number
                                sqc = sqcs, # quarantine scenarios
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
  my.files <- readRDS(paste0("output/MS_sim_runs_R0", R0, "_sqc", sqc, "_vc", vac, ".Rdata"))
  
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

# saveRDS(final_list_ts, "output/MS_monthly_infection_ts_10yrs.Rdata") 
# final_frame_box$pop[which(is.na(final_frame_box$pop))] <- 200000
# write.csv(final_frame_box, "output/MS_monthly_infection_boxplot_10yrs.csv", row.names=F) 
# write.csv(final_frame_box_burnout, "output/MS_monthly_infection_boxplot_burnout_10yrs.csv", row.names=F) 

#################################################################################################################################################
### BOXPLOT 1: Exposed and Infectious dogs per quarantine scenario
final_frame_box <- read.csv("output/MS_monthly_infection_boxplot_burnout_10yrs.csv") 
head(final_frame_box)

## colours
q4 <- sequential_hcl(4, "PuBuGn")
# plot(1:4, 1:4, col=q4, pch=16, cex=3)
mygray <- alpha("gray", 0.1)

## Set up theme
theme_set(theme_bw() +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.border=element_blank(),
                  axis.text=element_text(size=6),
                  axis.title.x = element_text(size = 9),
                  axis.title.y = element_text(size = 9),
                  # legend.title=element_blank(),
                  legend.text = element_text(size = 8),
                  text=element_text(size=9),
                  legend.key.size = unit(0.8, 'cm')))

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

ggsave("figs/InfectionBurninout_CompGrid.png", plot=Infection_CompGrid, width = 20, height = 9, units = "cm")

#################################################################################################################################################
### BOXPLOT 2: Incidence per 100,000 per quarantine scenario
## Draw plot
p_box <- ggplot(data = final_frame_box, aes(x=Quarantine, y=expD+infD/pop*10000, fill=Vaccination)) +
  geom_boxplot(outlier.size = 0.05, lwd=0.1, outlier.color=mygray) +
  # geom_boxplot(outlier.shape = NA) +
  facet_wrap(~R0, labeller = to_string, scales="free") +
  scale_fill_manual(name="Vaccination", values=c(q4), labels=c("0%", "25%", "50%", "75%")) + 
  ylab("Monthly incidence of infected dogs (E+I) / 10,000") + xlab("Quarantine scenario") +
  theme(strip.background = element_blank())

g1 <- p_box %+% dplyr::filter(final_frame_box, R0 == 1.3) + theme(legend.position = "none")
g2 <- p_box %+% dplyr::filter(final_frame_box, R0 != 1.3) + facet_wrap(~R0, nrow=2, scales="free")

Infection_CompGrid <- gridExtra::grid.arrange(g1, g2,
                                              layout_matrix = 
                                                matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                                         1, 1, 1, 2, 2, 2, 2, 2, 2),
                                                       byrow = TRUE, nrow = 2))

# ggsave("figs/Infection_CompGrid.png", plot=Infection_CompGrid, width = 20, height = 9, units = "cm")
ggsave("figs/IncidenceBurninout_CompGrid.png", plot=Infection_CompGrid, width = 20, height = 9, units = "cm")

#################################################################################################################################################
### BOXPLOT 3: Incidence per 100,000 per incursion scenario
final_frame_box <- read.csv("output/incs/MS_monthly_infection_boxplotinc_burnout_10yrs.csv") 
head(final_frame_box)

## Draw plot
tibble(final_frame_box)
final_frame_box <- tibble(final_frame_box) %>%
  mutate(R0 = as.factor(R0)) %>%
  mutate(Vaccination = as.factor(Vaccination)) %>%
  mutate(Incursion = as.factor(Incursion))

## Draw plot
p_box <- ggplot(data = final_frame_box, aes(x=Incursion, y=expD+infD/pop*10000, fill=Vaccination)) +
  geom_boxplot(outlier.size = 0.05, lwd=0.1, outlier.color=mygray) +
  # geom_boxplot(outlier.shape = NA) +
  facet_wrap(~R0, labeller = to_string, scales="free") +
  scale_fill_manual(name="Vaccination", values=c(q4), labels=c("0%", "25%", "50%", "75%")) + 
  ylab("Monthly incidence of infected dogs (E+I) / 10,000") + xlab("Incursion scenario") +
  theme(strip.background = element_blank())

g1 <- p_box %+% dplyr::filter(final_frame_box, R0 == 1.3) + theme(legend.position = "none")
g2 <- p_box %+% dplyr::filter(final_frame_box, R0 != 1.3) + facet_wrap(~R0, nrow=2, scales="free")

Incursion_CompGrid <- gridExtra::grid.arrange(g1, g2,
                                              layout_matrix = 
                                                matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                                         1, 1, 1, 2, 2, 2, 2, 2, 2),
                                                       byrow = TRUE, nrow = 2))

ggsave("figs/IncursionBurninout_CompGrid.png", plot=Incursion_CompGrid, width = 20, height = 9, units = "cm")



