## Quarantine model - thesis
## Code: weekly vaccination extrapolated for 5 years based on 2yr monthly data from Albay (2018 - 2019)

require(gamlss)
library("gamlss.dist")

vacc <- read.csv("~/Documents/Rabies_Warwick/Bicol_Philippines/BITERS/output/MonthlyVcBarangay_2014-2019.csv")
doses <- read.csv("~/Documents/Rabies_Warwick/Bicol_Philippines/BITERS/output/MonthlyVaccDosesBarangay_2014-2019.csv")
dogs <- read.csv("~/Documents/Rabies_Warwick/Bicol_Philippines/BITERS/output/MonthlyDogPopBarangay_2014-2019.csv")

### Get 2018 - 2019 & initial conditions
vacc<- vacc[,50:73]; dogs <- dogs[,50:73];doses <- doses[,50:73]
meancoverage <- colMeans(vacc)[1]
N <- colSums(dogs)[1]
V <- round(meancoverage*N, digits=0)
I <- round(49/56/0.05, digits=0) # assuming surveillance detected 5% of all cases in 13 mths=56 weeks
E <- round(4*I, digits=0)
V.comp <- rowSums(rmultinom(n=V, size=1, prob=c(50, 30, 20)))
V1 <- V.comp[1]
V2 <- V.comp[2]
V3 <- V.comp[3]
S <- round(as.numeric(N - I - E - V), digits=0)

### Calculate proportions weekly
props <- mapply(`/`, data.frame(doses), colSums(dogs))
propall <- colSums(props)
vaccweekly <- rep(propall, c(4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5))/4.3

# ### Look at distributions
# descdist(vaccweekly, discrete = F)
# # beta
# fit.beta <- fitdist(vaccweekly, distr = "beta", method = "mme"); plot(fit.beta)
# generated=rbeta(n=24, shape1=fit.beta$estimate[1], shape2=fit.beta$estimate[2])
# par(mfrow=c(2,1))
# hist(vaccweekly, breaks=seq(0,0.1,.001))
# hist(generated, breaks=seq(0,0.1,.001)) # but need zeros!

# # zero inflated beta
# library("gamlss.dist")
# require(gamlss)
# mod <- gamlss(vaccweekly~1,sigma.formula=~1, nu.formula=~1, family=BEZI) 
# generated=rBEZI(n=24, mu=fitted(mod,"mu")[1], sigma=fitted(mod,"sigma")[1], nu=fitted(mod,"nu")[1])
# par(mfrow=c(2,1))
# hist(vaccweekly, breaks=seq(0,0.1,.001))
# hist(generated, breaks=seq(0,0.1,.001))

### Draw weekly rates for 5 years from zero inflated Bet distribution
mod <- gamlss(vaccweekly~1,sigma.formula=~1, nu.formula=~1, family=BEZI) 
set.seed(1666318) # set seed
vaccweekly5yrs <- rBEZI(n=52*5, mu=fitted(mod,"mu")[1], sigma=fitted(mod,"sigma")[1], nu=fitted(mod,"nu")[1])

rm("vacc", "doses", "dogs", "mod", "props", "propall", "vaccweekly", 
   "meancoverage", "N", "V", "V.comp")


