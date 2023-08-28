## Quarantine model 
## Name: Isty Rysava
## Date: 15/02/2023
## Code: SEIRVQ deterministic model for rabies with carrying capacity 

SEIRVQ_rabies <- function(t, x, params) { 
  ## Input:
  ##  - t - time
  ##  - x - state variables
  ## Output: list of population sizes for each compartment of the model
  ## note: FD framework for the entire population
  
  # initialize
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  V <- x[5]
  Q <- x[6]
  
  N <- S + E + I + V + Q
  
  # calculate population sizes
  with(as.list(params), { 
    dS <- b*N - d*S - delta*N*S - beta*S*I/N - vc*S + wn*V 
    dE <- beta*S*I/N - d*E - delta*N*E - sigma*E 
    dI <- sigma*E - d*I - delta*N*I - gamma*I - q*I
    dR <- gamma*I + tau*Q + d*S + d*E + d*I + d*V + d*Q + delta*N*S + delta*N*E + delta*N*I + delta*N*V + delta*N*Q
    dV <- vc*S - d*V - delta*N*V - wn*V 
    dQ <- q*I - d*Q - delta*N*Q - tau*Q 
    
    res <- c(dS, dE, dI, dR, dV, dQ)
    list(res)
  }
  )
}

