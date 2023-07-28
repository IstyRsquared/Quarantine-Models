## Quarantine model 
## Name: Isty Rysava
## Date: 06/29/2022
## Code: Theoretical part - takes output from SEI(R)VQ stability analysis and draws plots

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

## Libraries 
library(harrypotter)
library(tidyverse)
library(gridExtra)
# library("ggpubr")

## Data
final_df <- read.csv("output/SEIRVQrabies_EndemicEquilibrium.csv")
unique(final_df$R0) # 11

## Prep col and labels
supp.labs <- paste("R0 =", unique(final_df$R0))
to_string <- as_labeller(c(`1` = supp.labs[1], `1.1` = supp.labs[2], `1.2` = supp.labs[3], `1.3` = supp.labs[4], `1.4` = supp.labs[5],
                           `1.5` = supp.labs[6], `1.6` = supp.labs[7], `1.7` = supp.labs[8], `1.8` = supp.labs[9], `1.9` = supp.labs[10], 
                           `2` = supp.labs[11]))

pal <- hp(n = 10, house = "Slytherin")
# plot(1:10, 1:10, col=pal, pch=16, cex=3)
colgreen <- pal[5]; collightgreen <- pal[2]

pal <- hp(n=8, option = "LunaLovegood")
# plot(1:8, 1:8, col=pal, pch=16, cex=3)
colpink <- pal[6]

colblue <- "dodgerblue4"
col <- colorRampPalette(c(colblue, colpink))(15)
# plot(1:15, 1:15, col=col, pch=16, cex=3)

# test <- filter(final_df, variable=="stability", R0=="1.5")
# unique(test$value)
# ggplot(test, aes( x=vc, y=q, fill = value)) +
#   geom_raster()
# y=unique(test$q)
# x=unique(test$vc)
# z=matrix(test$value, ncol=length(y), nrow=length(x))
# image(x=x, y=y, z=z, xlab="qP", ylab="vcP", col=col, main="Stability", cex.axis=.6, cex.main=.75, cex.lab=.7)

### SET UP
theme_set(theme_bw() +
            theme(plot.margin = unit(c(0.4, 0, 0.3, 0.1),
                                     "inches"),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.border=element_blank(),
                  axis.text=element_text(size=6),
                  axis.title.x = element_text(size = 9),
                  axis.title.y = element_text(size = 9),
                  legend.title=element_blank(),
                  legend.text = element_text(size = 8),
                  text=element_text(size=10),
                  legend.key.size = unit(0.8, 'cm')))

################# STABILITY PLOTS #################
p_stability <- ggplot(data = filter(final_df, variable=="stability"), aes(vc, q, fill = as.factor(value))) +
  geom_raster() +
  facet_wrap(~R0, labeller = to_string) +
  scale_fill_manual(values=c(colblue, colpink), name="", labels=c("stable", "unstable")) +
  xlab("Proportion of vaccinated dogs (vp)") + ylab("Proportion of quarantined dogs (qp)") + # labels need changing
  theme(strip.background = element_blank(),
        axis.title.y = element_text(vjust = +4),
        axis.title.x = element_text(vjust = -1))

stab_dat <- filter(final_df, variable=="stability")
g1 <- p_stability %+% dplyr::filter(stab_dat, R0 == 1.3) + theme(legend.position = "none")
g2 <- p_stability %+% dplyr::filter(stab_dat, R0 != 1.3) + facet_wrap(~R0, nrow=2)

stability_grid <- gridExtra::grid.arrange(g1, g2,
                        layout_matrix = 
                          matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                   1, 1, 1, 2, 2, 2, 2, 2, 2),
                                 byrow = TRUE, nrow = 2))

# stability_grid <- gridExtra::grid.arrange(g1, g2,
#                                           layout_matrix = 
#                                             matrix(c(1, 1, 1, NA, NA, NA, NA, NA, NA,
#                                                      1, 1, 1, 2, 2, 2, 2, 2, 2,
#                                                      1, 1, 1, 2, 2, 2, 2, 2, 2,
#                                                      1, 1, 1, NA, NA, NA, NA, NA, NA),
#                                                    byrow = TRUE, nrow = 4))

################# INFECTION PLOTS #################
### A) Infected dogs
inf_dat <- filter(final_df, variable=="infected")
inf.dogs <- tibble(inf_dat)%>%
  mutate(value=ifelse(value<0, NA, value)) # for unstable NA, rest by number!
inf.dogs$inc_10t <- inf.dogs$value/filter(final_df, variable=="pop")$value*10000

p_popinf <- ggplot(data = inf.dogs, aes(vc, q, fill = inc_10t)) +
  geom_raster() +
  facet_wrap(~ R0, labeller = to_string) +
  xlab("Proportion of vaccinated dogs (vp)") + ylab("Proportion of quarantined dogs (qp)") +
  scale_fill_gradientn(colours=col, na.value = "gray90") +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(vjust = +4),
        axis.title.x = element_text(vjust = -1))

g1 <- p_popinf %+% dplyr::filter(inf.dogs, R0 == 1.3) + theme(legend.position = "none")
g2 <- p_popinf %+% dplyr::filter(inf.dogs, R0 != 1.3) + facet_wrap(~R0, nrow=2)

infection_grid <- gridExtra::grid.arrange(g1, g2,
                                          layout_matrix = 
                                            matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                                     1, 1, 1, 2, 2, 2, 2, 2, 2),
                                                   byrow = TRUE, nrow = 2))
### B) Exposed dogs
exp_dat <- filter(final_df, variable=="exposed")
exp.dogs <- tibble(exp_dat)%>%
  mutate(value=ifelse(value<0, NA, value)) # for unstable NA, rest by number!
exp.dogs$inc_10t <- exp.dogs$value/filter(final_df, variable=="pop")$value*10000

p_popexp<- ggplot(data = exp.dogs, aes(vc, q, fill = inc_10t)) +
  geom_raster() +
  facet_wrap(~ R0, labeller = to_string) +
  xlab("Proportion of vaccinated dogs (vp)") + ylab("Proportion of quarantined dogs (qp)") +
  scale_fill_gradientn(colours=col, na.value = "gray90") +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(vjust = +4),
        axis.title.x = element_text(vjust = -1))

g1 <- p_popexp %+% dplyr::filter(exp.dogs, R0 == 1) + theme(legend.position = "none")
g2 <- p_popexp %+% dplyr::filter(exp.dogs, R0 != 1) + facet_wrap(~R0, nrow=2)

exposed_grid <- gridExtra::grid.arrange(g1, g2,
                                          layout_matrix = 
                                            matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2,
                                                     1, 1, 1, 2, 2, 2, 2, 2, 2),
                                                   byrow = TRUE, nrow = 2))
## Bind
# Stability_Infection_Grid <- ggpubr::ggarrange(stability_grid, infection_grid,
#                                               ncol=1, labels=c("Stability","Infection"),
#                                               hjust = -0.1,
#                                               font.label = list(size = 5))

Stability_InfectedPop_Grid <- ggpubr::ggarrange(stability_grid, infection_grid, exposed_grid,
                                            ncol=1, 
                                            labels=c("A) Stability","B) Infected incidence/10,000", "C) Exposed incidence/10,000"),
                                            hjust = -0.1,
                                            font.label = list(size = 10))

ggsave("figs/Stability_InfectedPop_Grid.png", width = 22, height = 28, units = "cm")


# ggsave("figs/Stability_Infection_Grid.png", width = 20, height = 18, units = "cm")

