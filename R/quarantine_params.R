## Computational part - SEI(R)VQ tauleap model
## Code: parameters to test quarantine scenarios

### SCENARIO 1 - no quarantine 
## q1, range 1, q2, range 2 = 0
sqc1 <- c(b=0.38/52, d=0.28/52,   
             R0=1.3, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7,vw=1/52, vv=0.4, i=1.5,
             q1=0, range1=0, q2=0, range2=0, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 2 - quarantine of incident-reported dogs
## q1=2.5471, range1=4.3995, q2=2.5, range2=4.3 
sqc2 <- c(b=0.38/52, d=0.28/52,   
             R0=1.2, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=2.5, range2=4.3, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)

### SCENARIO 3 - quarantine of incident-reported dogs + follow up contact tracing
## q1=2.5471, range 1=4.3995, q2=2.8, range2=1.5
sqc3 <- c(b=0.38/52, d=0.28/52,   
             R0=1.2, size=1.33, Rprob=0, sigma= 1/22.3*7, gamma=1/7*7, vw=1/52, vv=0.4, i=1.5,
             q1=2.5471, range1=4.3995, q2=3, range2=1.3, eta=0.7, qr=1/14*7, 
             bh=0.019/52, dh=0.0052/52, hrr=0.2, meanpat=1.259512e-04, sdpat=5.479777e-05,
             R0h=0.37, sizeh=0.56, Rprobh=0.56/(0.37+0.56), sigmah=1/40*7, gammah=1/7*7, vs=1/8*7, mue=0.8, lambda=1/5*7)
params_list <- list(sqc1=sqc1, sqc2=sqc2, sqc3=sqc3)
                    
# # Test R0 distributions
# k=1.33
# R0s=seq(1, 2, 0.1) 
# 
# pdf("figs/R0parmdist_exploreMean.pdf", width=4, height=8)
# par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
# for(R0 in R0s){
#     p <- k/(R0+k)
#     dist.tmp <- rnbinom(1000, size=k, prob=p)
#     hist(dist.tmp, breaks=seq(-1, max(dist.tmp)+1, by=1),
#          main=paste("R0:", R0, "p:", p), cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
# }
# dev.off()
# 
# pdf("figs/R0parmdist_exploreMu.pdf", width=4, height=8)
# par(mfrow=c(6,3), mar=c(2.5,2,1.5,1))
# for(R0 in R0s){
#   dist.tmp <- rnbinom(1000, mu=R0, size=k)
#   hist(dist.tmp, breaks=seq(-1, max(dist.tmp)+1, by=1),
#        main=paste("R0:", R0), cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
# }
# dev.off()
# 
# R0=0.3; k=0.56
# p <- k/(R0+k)
# dist.tmp1 <- rnbinom(1000, size=k, prob=p)
# hist(dist.tmp, breaks=seq(-1, max(dist.tmp1)+1, by=1),
#  main=paste("R0:", R0, "p:", p), cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
# 
# dist.tmp2 <- rnbinom(1000, mu=R0, size=k)
# hist(dist.tmp, breaks=seq(-1, max(dist.tmp2)+1, by=1), 
#      main=paste("R0:", R0), cex.main=0.5, cex.lab=0.5, cex.axis=0.5)
# table(dist.tmp1)
# table(dist.tmp2)





