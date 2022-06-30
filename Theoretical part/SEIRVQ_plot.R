## Quarantine model 
## Name: Isty Rysava
## Date: 06/29/2022
## Code: Theoretical part - takes output from SEI(R)VQ stability analysis and draws plots

rm(list=ls())
setwd("~/Documents/Rabies_Warwick/Quarantine-models")

## Libraries 
library(harrypotter)

## Data
final_df <- read.csv("output/SEIRVQrabies_EndemicEquilibrium.csv")

## Name: Micah Fletcher (adapted by Isty Rysava)
## Date: 21/06/2022
## Code: Makes faceted raster plots for Anthrax_IBM output data

rm(list=ls())
setwd("~/Documents/Anthrax-Virulence-main")

# LIBRARIES
library("tidyverse")
library("ggpubr")
library("harrypotter")

# DATA
final_df <- read.csv("output/results_summary.csv")

df <- tibble(final_df) %>%
  mutate(surv = ifelse(is.na(surv), "Full", "Half"),
         Kvar = ifelse(Kvar==TRUE, "Fluctuate", "Constant"),
         omega = factor(omega), d = factor(Basemorts))

### SET UP
theme_set(theme_bw() +
            theme(panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  panel.border=element_blank(),
                  axis.text=element_text(size=5),
                  axis.title.x = element_text(size = 5),
                  axis.title.y = element_text(size = 5),
                  legend.title=element_blank(),
                  legend.text = element_text(size = 4),
                  text=element_text(size=5),
                  legend.key.size = unit(0.3, 'cm')))

################# PREP PLOTS 1 #################
p_toxin <- ggplot(data = filter(df, variable=="toxin", surv=="Full"), aes(omega, d, fill = mean)) +
  geom_raster() +
  facet_grid(Kvar ~ K0) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="yellow", 
                       midpoint=600, limits=range(filter(df, variable=="toxin", surv=="Full")$mean)) +
  theme(strip.background = element_blank())

p_prt <- ggplot(data = filter(df, variable=="pathogen", surv=="Full"), aes(omega, d, fill = mean)) +
  geom_raster() +
  facet_grid(Kvar ~ K0) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="yellow", 
                       midpoint=3E7, limits=range(filter(df, variable=="pathogen", surv=="Full")$mean)) +
  theme(strip.background = element_blank())

p_infD <- ggplot(data = filter(df, variable=="infD", surv=="Full"), aes(omega, d, fill = mean)) +
  geom_raster() +
  facet_grid(Kvar ~ K0) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="yellow", 
                       midpoint=70, limits=range(filter(df, variable=="infD", surv=="Full")$mean)) +
  theme(strip.background = element_blank())

p_strat <- ggplot(data = filter(df, variable=="strategy", surv=="Full"), aes(omega, d, fill = mean)) +
  geom_raster() +
  facet_grid(Kvar ~ K0) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="yellow", 
                       midpoint=0.5, limits=range(filter(df, variable=="strategy", surv=="Full")$mean)) +
  theme(strip.background = element_blank())

### DRAW PLOTS
(summary_plot_survFull <- ggpubr::ggarrange(p_toxin, p_prt, p_infD, p_strat,
                                            ncol=1, labels=c("Toxin load","Pathogen count","Infection duration", "Strategy"),
                                            hjust = -0.1,
                                            font.label = list(size = 5)))

# ggsave(faceted_raster_plot, device="png")
ggsave("summary_plot_survFull.png", width = 10, height = 17, units = "cm")
# ggsave("faceted_raster_plotA4.png", width = 21, height = 29, units = "cm")