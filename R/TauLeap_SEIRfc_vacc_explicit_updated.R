## Quarantine model 
## Name: Isty Rysava
## Date: 06/29/2022
## Code: Computational part - SEI(R)VQ tauleap model function
##  features: biased vaccination, incursions, explicit biting behaviour, quarantine 

library(compoisson)
library(truncnorm)

# N.B. go back to the original format for HK model!
SEIR.tauleap <- function(init, pars, end.time, tau){
  init2 <- init
  
  Equations <- function(params, init, end.time, tau) {
    with(as.list(c(pars, init)), {
      x <- init
      tau <- tau
      rate <- rep(0, 33) 
      counts <- rep(0, length(init)) 
      change <- matrix(0, nrow=33, ncol=length(x))
      
      ## previous state
      S.old=x[1]
      V1.old=x[2]
      V2.old=x[3]
      V3.old=x[4]
      E.old=x[5]
      I.old=x[6]
      Qs.old=x[7]
      Qer.old=x[8] # quarantined exposed dogs that will be REMOVED (die)
      Qeb.old=x[9] # quarantined exposed dogs that will go BACK into community (is this leaky quarantine)
      Qi.old=x[10]
      Rill.old=x[11]
      
      Sh.old=x[12]
      Eh.old=x[13]
      Ih.old=x[14]
      Rh.old=x[15]
      Vhs.old=x[16]
      Vhe.old=x[17]
      
      N = S.old + V1.old + V2.old + V3.old + E.old + I.old  
      Nh = Sh.old + Eh.old + Ih.old + Vhs.old + Vhe.old
      
      ## Set parameters:
      b <- params["b"] # birth rate dogs
      d <- params["d"] # death rate dogs
      sigma <- params["sigma"] # 1/incubation period dogs
      gamma <- params["gamma"] # recovery/removal rate dogs; 1/infectious period
      vc <- params["vc"] # vaccination dogs
      vw <- params["vw"] # waning immunity dogs
      vv <- params["vv"] # vaccination bias
      qr <- params["qr"] # rate of quarantine; 1/duration of quarantine
      i <- params["i"] # incursion rate
      
      bh <- params["bh"] # birth rate humans
      dh <- params["dh"] #  death rate human
      sigmah <- params["sigmah"] # 1/incubation period humans
      gammah <- params["gammah"] # 1/infectious period humans
      vs <- params["vs"] # vacciation rate humans; 1/duration it takes to achieve immunity (assumed 2 doses would be enough)
      mue <- params["mue"] # % of exposed humans presenting at clinics
      lambda <- params["lambda"] # 1/period after which patients are removed from the V compartment (bc dog quarantine is based of this number they cannot stay there forever)
      hrr <- params["hrr"] # risk of developing rabies in humans upon rabies exposure <- 0.2
      epsilon <- rtruncnorm(1, a=0, b=1, mean=params["meanpat"], sd=params["sdpat"]) # proportion of the total population showing up at clinics
      
      ## Calculate rates & define events
      rate[1] <- b*N # birth.rate
      change[1, ] <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[2] <- d*S.old # deathS.rate
      change[2, ] <- c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[3] <- d*V1.old # deathV1.rate
      change[3, ] <- c(0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[4] <- d*V2.old # deathV2.rate
      change[4, ] <- c(0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[5] <- d*V3.old # deathV3.rate
      change[5, ] <- c(0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[6] <- d*E.old # deathE.rate
      change[6, ] <- c(0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      
      # maybe here I need exposed back to sus 50% chance; -log(1-0.5)/2.5 weeks
      rate[7] <- sigma*E.old # sigma.rate
      change[7, ] <- c(0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[8] <- gamma*I.old # gamma.rate
      change[8, ] <- c(0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
      
      rate[9] <- sigma*Qer.old # sigma.q
      change[9, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0)
      rate[10] <- gamma*Qi.old # gamma.q
      change[10, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
      
      rate[11] <- qr*Qs.old # qr.S 
      change[11, ] <- c(1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[12] <- qr*Qer.old # qr.Er 
      change[12, ] <- c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0)
      rate[13] <- qr*Qeb.old # qr.Eb 
      change[13, ] <- c(0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[14] <- qr*Qi.old # qr.I 
      change[14, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0)
      
      rate[15] <- vc*(1-vv)*S.old # vcS.rate
      change[15, ] <- c(-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[16] <- vc*V1.old + vc*vv*S.old*(1/3) # vcV1.rate
      change[16, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[17] <- vc*V2.old + vc*vv*S.old*(1/3) # vcV2.rate
      change[17, ] <- c(0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[18] <- vc*V3.old + vc*vv*S.old*(1/3) # vcV3.rate 
      change[18, ] <- c(0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[19] <- vw*V1.old # vw.V1 
      change[19, ] <- c(0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[20] <- vw*V2.old # vw.V2 
      change[20, ] <- c(0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      rate[21] <- vw*V3.old # vw.V3 
      change[21, ] <- c(1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      
      ## HUMANS
      rate[22] <- bh*Nh # birth.h
      change[22, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
      rate[23] <- dh*Sh.old # death.Sh
      change[23, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0)
      rate[24] <- dh*Eh.old # death.Eh
      change[24, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0)
      rate[25] <- dh*Vhs.old # death.Vhs
      change[25, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0)
      rate[26] <- dh*Vhe.old # death.Vhe
      change[26, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1)
      
      rate[27] <- sigmah*Eh.old # sigma.h 
      change[27, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0)
      rate[28] <- gammah*Ih.old # gamma.h
      change[28, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0)
      
      rate[29] <- Eh.old*mue*vs # vcVhe (vacc E) # mue = % of exposed coming to clinics 
      change[29, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1)
      rate[30] <- Sh.old*vs*epsilon # vcVhs (vacc S)
      change[30, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0)
      
      rate[31] <- Eh.old*(1-hrr)*(sigmah*2) # exposed to sus
      change[31, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0)
      rate[32] <- Vhe.old*lambda # Vhe to leave
      change[32, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1) # move to sus
      rate[33] <- Vhs.old*lambda # Vhs to leave
      change[33, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0) # move to sus
      
      ## Tau leap
      count.change <- change
      count.change[count.change<0] <- 0
      for (j in 1:length(rate)) {
        num <- rpois(1, rate[j] * tau)
        num.min <- min(num, init[which(change[j, ] < 0)])
        init <- init + change[j, ] * num.min
        counts <- counts + count.change[j, ] * num.min
      }
      
      names(init) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh", 
                       "Ih", "Rh", "Vhs", "Vhe")
      names(counts) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh", 
                         "Ih", "Rh", "Vhs", "Vhe")
      
      ## currents state
      y <- init
      S=y[1]
      V1=y[2]
      V2=y[3]
      V3=y[4]
      E=y[5]
      I=y[6]
      Qs=y[7]
      Qer=y[8]
      Qeb=y[9]
      Qi=y[10]
      Rill=y[11]
      
      Sh=y[12]
      Eh=y[13]
      Ih=y[14]
      Rh=y[15]
      Vhs=y[16]
      Vhe=y[17]
      
      ## Explicitly defined parameters: 
      ## 1) R0 dogs and humans
      # dogs R0
      if(I>0){
        # print(paste0("I", I))
        Rnull.temp <- rnbinom(I, size=params["size"], prob=params["Rprob"])
        # Rnull.temp <- rnbinom(I, mu=params["R0"], size=1/params["size"])
        # print(paste0("new cases", sum(Rnull.temp)))
        
      }else{
        Rnull.temp <- 0
      }
      
      if(sum(Rnull.temp)!=0){ 
        new.bitestemp <- sum(Rnull.temp)*tau*params["gamma"] 
        new.bites <- new.bitestemp
        # print(paste0("new bites temp:", new.bitestemp))
        # if(runif(1)<new.bitestemp){
        #   new.bites <- floor(sum(Rnull.temp)*tau*params["gamma"])+1
        # }else{
        #   new.bites <- floor(sum(Rnull.temp)*tau*params["gamma"])
        # }
        if((S+V1+V2+V3+E+I)>new.bites){
          pop.per.compartment <- c(S, V1, V2, V3, E, I)
          comp.bitten <- rowSums(rmultinom(n=new.bites, size=1, prob = pop.per.compartment))
          S <- S - comp.bitten[1]
          V1 <- V1
          V2 <- V2
          V3 <- V3
          E <- E + comp.bitten[1]
          I <- I 
          counts["E"] <- comp.bitten[1]
        }else{
          print("this is crazy dogs")
        }
      }
      
      # humans R0: 
      if(I>0){
        Rnullh.temp <- rnbinom(I, size=params["sizeh"], prob=params["Rprobh"])
      }else{
        Rnullh.temp <- 0
      }
      
      if (sum(Rnullh.temp)!=0){ 
        new.bitestemp <- sum(Rnullh.temp)*tau*params["gamma"] 
        new.bites <- new.bitestemp
        # if(runif(1)<new.bitestemp){
        #   new.bites <- floor(sum(Rnullh.temp)*tau*params["gamma"])+1
        # }else{
        #   new.bites <- floor(sum(Rnullh.temp)*tau*params["gamma"])
        # }
        if (Sh>new.bites){
          Sh <- Sh - new.bites
          Eh <- Eh + new.bites
        }else{
          print("this is crazy humans")
        }
      }
      
      ## 2) quarantine
      ## positive quarantine rate
      if(sum(params["q2"], params["range2"])!=0 & Vhe>0){
        qdp <- rcom(Vhe, lambda=params["q2"], nu=params["range2"]) 
      }else{
        qdp <- 0
      }
      
      if (sum(qdp)!=0){
        new.qarnt <- sum(qdp)*tau - floor(sum(qdp)*tau)
        if(runif(1)<new.qarnt){
          q.dogs <- floor(sum(qdp)*tau)+1
        }else{
          q.dogs <- floor(sum(qdp)*tau)
        }
        if(q.dogs>0){
          # remove biting Is (culprit) - move to Qi
          sumq.dogs <- q.dogs - length(which(qdp>0))
          Qi <- Qi + min(I, length(which(qdp>0))) # qr
          I <- max(0, I - length(which(qdp>0))) # inf
          
          # distribute remaining
          case.dogs <- round(sumq.dogs * params["eta"], digits=0) # eta=0.8 (proportion of dogs identified through follow-up as suspect)
          noncase.dogs <- sumq.dogs-case.dogs
          
          # cases
          disease.periods <- c(3.1, 22.3) # infectious period, incubation period
          I.tmp <- rowSums(rmultinom(n=case.dogs, size=1, prob = disease.periods))[1] # distribute proportionally to the duration of each period
          E.tmp <- case.dogs-I.tmp
          Er.tmp <- ifelse(E.tmp>0, length(which(runif(E.tmp, 1,22.3)<(7/params["qr"]))), 0) # allow for leaky quarantine; exposed that will become rabid during 14-day observation
          # if(E.tmp>0){
          #   Er.tmp <- length(which(runif(E.tmp, 1,28)<(7/params["qr"])))
          # }else{
          #   Er.tmp <- 0
          # }
          Eb.tmp <- E.tmp-Er.tmp # leaky quarantine
          
          Qi <- Qi + min(I, I.tmp)
          I <- max(0, I - I.tmp)
          
          Qer <- Qer + min(E, Er.tmp)
          E <- max(0, E - Er.tmp)
          Qeb <- Qeb + min(E, Eb.tmp)
          E <- max(0, E - Eb.tmp)
          
          # non-cases
          pop.per.compartment <- c(S, V1, V2, V3)
          comp.q <- rowSums(rmultinom(n=noncase.dogs, size=1, prob = pop.per.compartment))
          S <- S - comp.q[1]
          Qs <- Qs + comp.q[1]
          V1 <- V1 
          V2 <- V2 
          V3 <- V3 
        }
      }
      
      ## false quarantine rate
      if(sum(params["q1"], params["range1"])!=0 & Vhs>0){
        qdf <- rcom(Vhs, lambda=params["q1"], nu=params["range1"]) 
      }else{
        qdf <- 0
      }
      
      if (sum(qdf)!=0){
        new.qarnt <- sum(qdf)*tau - floor(sum(qdf)*tau)
        if(runif(1)<new.qarnt){
          q.dogs <- floor(sum(qdf)*tau)+1
        }else{
          q.dogs <- floor(sum(qdf)*tau)
        }
        if(q.dogs>0){
          # all are non-cases
          pop.per.compartment <- c(S, V1, V2, V3, E)
          comp.q <- rowSums(rmultinom(n=q.dogs, size=1, prob = pop.per.compartment))
          S <- S - comp.q[1]
          Qs <- Qs + comp.q[1]
          V1 <- V1 
          V2 <- V2 
          V3 <- V3 
          
          if(comp.q[5]>0){
            Er.tmp <- length(which(runif(comp.q[5], 1, 22.3)<(7/params["qr"]))) # allow for leaky quarantine; exposed that will become rabid during 14-day observation
            Eb.tmp <- comp.q[5]-Er.tmp # these will go back to I
            E <- E - comp.q[5]
            Qer <- Qer + Er.tmp # removed
            Qeb <- Qeb + Eb.tmp # back
          }
        }
      }
      
      ## #) Incursions 
      incs <- rpois(1, i)
      I <- I+incs
      counts["I"] <- counts["I"]+incs
      
      init <- c(S=S, V1=V1, V2=V2, V3=V3, E=E, I=I, Qs=Qs, Qer=Qer, Qeb=Qeb, Qi=Qi, Rill=Rill,
                Sh=Sh, Eh=Eh, Ih=Ih, Rh=Rh, Vhs=Vhs, Vhe=Vhe)
      names(init) <- c("S", "V1", "V2", "V3", "E", "I", "Qs", "Qer", "Qeb", "Qi", "Rill", "Sh", "Eh", 
                       "Ih", "Rh", "Vhs", "Vhe")
      
      return(list(init, counts))
    })
  }
  
  S <- V1 <- V2 <- V3 <- E <- I <- Qs <- Qer <- Qeb <- Qi <- Rill <- Sh <- Eh <- Ih <- Rh <- Vhs <- Vhe <- double()
  mycount <- c()
  t <- 0
  time <- seq(0, end.time, by = tau)
  for (t in 1:length(time)) {
    params <- pars
    tmp <- Equations(params, init, end.time, tau=tau)
    S <- c(S, init['S'])
    V1 <- c(V1, init['V1'])
    V2 <- c(V2, init['V2'])
    V3 <- c(V3, init['V3'])
    E <- c(E, init['E'])
    I <- c(I, init['I'])
    Qs <- c(Qs, init['Qs'])
    Qer <- c(Qer, init['Qer'])
    Qeb <- c(Qeb, init['Qeb'])
    Qi <- c(Qi, init['Qi'])
    Rill <- c(Rill, init['Rill'])
    Sh <- c(Sh, init['Sh'])
    Eh <- c(Eh, init['Eh'])
    Ih <- c(Ih, init['Ih'])
    Rh <- c(Rh, init['Rh'])
    Vhs <- c(Vhs, init['Vhs'])
    Vhe <- c(Vhe, init['Vhe'])
    init <- tmp[[1]]
    mycount <- c(mycount, tmp[[2]])
  }
  return(list(pars = pars,
              init = init2,
              time = time,
              counts = matrix(mycount, ncol=length(init), byrow=T),
              results = data.frame(time, S, V1, V2, V3, E, I, Qs, Qer, Qeb, Qi, Rill, 
                                   Sh, Eh, Ih,  Rh, Vhs, Vhe)))
}

