# Hong KONG model development
## Code: weekly vaccination extrapolated for 5 years based on 2yr monthly data from Albay

vacc <- read.csv("~/Dropbox/HongKong/data/AlbayVaccPopulationProp_bymonthbybarangay_Jan2015-Dec2016.csv")

### VACC 1 - simple repetition
vaccall <- colSums(vacc[,2:25])
vaccweekly <- rep(vaccall, c(4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly <- rep(vaccweekly,3)[1:(5*52)]
vacc1 <- vaccweekly

### VACC 2 - reshuffle months for each year
vaccall1 <- colSums(vacc[,2:13])
vaccall2 <- colSums(vacc[,14:25])
set.seed(1666318) # set seed
vaccall3 <- sample(vaccall1)
vaccall4 <- sample(vaccall2)
vaccall5 <- sample(vaccall1)

vaccweekly1 <- rep(vaccall1, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly2 <- rep(vaccall2, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly3 <- rep(vaccall3, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly4 <- rep(vaccall4, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly5 <- rep(vaccall5, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3

vaccweekly <- c(vaccweekly1, vaccweekly2, vaccweekly3, vaccweekly4, vaccweekly5)
vacc2 <- vaccweekly

### VACC 3 - resufle weeks each month
vaccall1 <- colSums(vacc[,2:13])
vaccall2 <- colSums(vacc[,14:25])

vaccweekly1 <- rep(vaccall1, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
vaccweekly2 <- rep(vaccall2, c(4,4,5,4,4,5,4,4,5,4,4,5))/4.3
set.seed(1666318) # set seed
vaccweekly3 <- sample(vaccweekly1)
vaccweekly4 <- sample(vaccweekly2)
vaccweekly5 <- sample(vaccweekly1)

vaccweekly <- c(vaccweekly1, vaccweekly2, vaccweekly3, vaccweekly4, vaccweekly5)
vacc3 <- vaccweekly

rm("vaccall", "vaccall1", "vaccall2", "vaccall3", "vaccall4", "vaccall5",
   "vaccweekly", "vaccweekly1", "vaccweekly2", "vaccweekly3", "vaccweekly4", "vaccweekly5")



