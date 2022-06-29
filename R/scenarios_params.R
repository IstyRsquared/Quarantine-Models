## Hong KONG model development
## Code: parameters to test quarantine scenarios

library(compoisson)

### SCENARIO 1 - no quarnatine or follow up
## q1, range 1, q2, range 2 = 0, tsn = 0, tsp = 1
params1 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, vw=1/52, vv=0.4, i=1.5,
             q1=0, range1=0, q2=0, range2=0, eta=0.8, qr=1/14*7, tsn=0, tsp=1,
             bh=0.019/52, dh=0.0052/52, hrr=0.2,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/10*7, mue=0.7, lambda=1/5*7)

## NOTE: vs assumed to be 10 day as the time from the bite incidet to the second/third PEP day (but could be longer depending on circumstances/settings)

### SCENARIO 2 - quarantine of incident-reported dogs
## q1=2.5471, range1=4.3995, q2=2.5, range2=4.3 , tsn = 0.85, tsp = 1
params2 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=2.5, range2=4.3, eta=0.8, qr=1/14*7, tsn=1, tsp=1,
             bh=0.019/52, dh=0.0052/52, hrr=0.2,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/10*7, mue=0.7, lambda=1/5*7)

### SCENARIO 3 - quarantine of incident-reported dogs + follow up
## q1=2.5471, range 1=4.3995, q2=2.8, range2=1.5, tsn = 0.85, tsp = 1
params3 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=3, range2=1.3, eta=0.8, qr=1/14*7, tsn=1, tsp=1,
             bh=0.019/52, dh=0.0052/52, hrr=0.2,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/10*7, mue=0.7, lambda=1/5*7)

### SCENARIO 4 - quarantine of incident-reported dogs + follow up + patient interviews 
## (bite & vacc histories)
## same as scenario 3, but calculate PEP as 50% of SUS get 1 dose and 50% 2 doses
# hist(rcom(1000, lambda=3, nu=1.4), breaks=seq(-1,8,1))
# hist(rcom(1000, lambda=2.5, nu=4.3), breaks=seq(-1,8,1))
# hist(rcom(1000, lambda=3, nu=1.3), breaks=seq(-1,8,1))



