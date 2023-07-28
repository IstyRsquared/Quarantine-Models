## Computational part - SEI(R)VQ tauleap model
## Code: parameters to test incursion scenarios

### SCENARIO 1 - no incursions
sqc1 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7,vw=1/52, vv=0.4, i=0,
             q1=0, range1=0, q2=0, range2=0, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 2 - incursion rate = 1.5
sqc2 <- c(b=0.38/52, d=0.28/52,   
          R0=1.3, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7,vw=1/52, vv=0.4, i=1.5,
          q1=0, range1=0, q2=0, range2=0, eta=0.7, qr=1/14*7, 
          bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
          R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 3 - incursion rate = 3
sqc3 <- c(b=0.38/52, d=0.28/52,   
          R0=1.3, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7,vw=1/52, vv=0.4, i=3,
          q1=0, range1=0, q2=0, range2=0, eta=0.7, qr=1/14*7, 
          bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
          R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)
params_list <- list(sqc1=sqc1, sqc2=sqc2, sqc3=sqc3)
                    
