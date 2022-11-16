## Quarantine model 
## Name: Isty Rysava
## Date: 06/29/2022
## Code: Theoretical part - takes output from SEI(R)VQ stability analysis and draws plots

rm(list=ls())
setwd("C:/Users/tui9/Documents/Practice code/Quarantine-Models")

## Libraries 
library(harrypotter)
library("tidyverse")
# library("ggpubr")

## Data
final_df <- read.csv("output/SEIRVQrabies_EndemicEquilibrium.csv")
unique(final_df$R0) # 11

## test
pal <- hp(n = 10, house = "Slytherin")
plot(1:10, 1:10, col=pal, pch=16, cex=3)
colgreen <- pal[5]
collightgreen <- pal[2]

pal <- hp(n=8, option = "LunaLovegood")
plot(1:8, 1:8, col=pal, pch=16, cex=3)
colpink <- pal[6]

col <- colorRampPalette(c(colgreen, colpink))(30) 
col.more <- colorRampPalette(c(collightgreen, colpink))(60)
plot(1:30, 1:30, col=col, pch=16, cex=3)
plot(1:60, 1:60, col=col.more, pch=16, cex=3)

test <- filter(final_df, variable=="stability", R0=="1.5")
unique(test$value)
ggplot(test, aes( x=vc, y=q, fill = value)) +
  geom_raster()

y=unique(test$q)
x=unique(test$vc) # the order doesn't match!
z=matrix(test$value, ncol=length(y), nrow=length(x))
image(x=x, y=y, z=z, xlab="qP", ylab="vcP", col=col.more, main="Stability", cex.axis=.6, cex.main=.75, cex.lab=.7)


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
p_stability <- ggplot(data = filter(final_df, variable=="stability"), aes(vc, q, fill = as.factor(value))) +
  geom_raster() +
  facet_grid(~R0) +
  scale_fill_manual(values=c(colgreen, colpink), name="") +
  xlab("Vaccination rate (vc)") + ylab("Quarantine rate (q)") +
  theme(strip.background = element_blank())

p_popinf <- ggplot(data = filter(final_df, variable=="infection"), aes(vc, q, fill = value)) +
  geom_raster() +
  facet_grid(~ R0) +
  scale_fill_gradient2(low="dodgerblue4", mid="white", high="yellow", 
                       midpoint=0, limits=range(filter(final_df, variable=="infection")$value)) +
  xlab("Vaccination rate (vc)") + ylab("Quarantine rate (q)") +
  theme(strip.background = element_blank())

### DRAW PLOTS
(summary_plot_survFull <- ggpubr::ggarrange(p_toxin, p_prt, p_infD, p_strat,
                                            ncol=1, labels=c("Stability","Infected population"),
                                            hjust = -0.1,
                                            font.label = list(size = 5)))

# ggsave(faceted_raster_plot, device="png")
ggsave("summary_plot_survFull.png", width = 10, height = 17, units = "cm")
# ggsave("faceted_raster_plotA4.png", width = 21, height = 29, units = "cm")