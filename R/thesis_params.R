## Quarantine model - thesis
## Code: parameters to test quarantine scenarios

### SCENARIO 1 - no quarnatine 
## q1, range 1, q2, range 2 = 0
params1 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, gammaR=1/3.1*7,vw=1/52, vv=0.4, i=1.5,
             q1=0, range1=0, q2=0, range2=0, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 2 - quarantine of incident-reported dogs
## q1=2.5471, range1=4.3995, q2=2.5, range2=4.3 
params2 <- c(b=0.38/52, d=0.28/52,   
             R0=1.2, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, gammaR=1/3.1*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=2.5, range2=4.3, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 3 - quarantine of incident-reported dogs + follow up contact tracing
## q1=2.5471, range 1=4.3995, q2=2.8, range2=1.5
params3 <- c(b=0.38/52, d=0.28/52,   
             R0=1.2, size=1.33, sigma= 1/22.3*7, gamma=1/7*7, gammaR=1/3.1*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=3, range2=1.3, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)




