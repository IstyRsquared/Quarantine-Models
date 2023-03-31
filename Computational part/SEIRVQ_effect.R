## Quarantine model 
## Name: Isty Rysava
## Date: 31/03/2023
## Code: Statistical analysis of vacination and quarantine effect on the number of cases



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
