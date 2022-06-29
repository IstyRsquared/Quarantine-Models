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
