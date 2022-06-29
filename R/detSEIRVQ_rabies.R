## Quarantine model 
## Name: Isty Rysava
## Date: 15/02/2021
## Code I: SEI(R)VQ model for rabies: (A) with carrying capacity, (B) with carrying capacity & reactive quarantine
## Code II: Newton function for equilibria of vector field

## ## Code I (A)
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

SEIVQ_rabies <- function(t, x, params) { 
  ## Input:
  ##  - t - time
  ##  - x - state variables
  ## Output: list of population sizes for each compartment of the model
  ## note: FD framework for the entire population
  
  # initialize
  S <- x[1]
  E <- x[2]
  I <- x[3]
  V <- x[4]
  Q <- x[5]
  
  N <- S + E + I + V + Q
  
  # calculate population sizes
  with(as.list(params), { 
    dS <- b*N - d*S - delta*N*S - beta*S*I/N - vc*S + wn*V 
    dE <- beta*S*I/N - d*E - delta*N*E - sigma*E 
    dI <- sigma*E - d*I - delta*N*I - gamma*I - q*I
    dV <- vc*S - d*V - delta*N*V - wn*V 
    dQ <- q*I - d*Q - delta*N*Q - tau*Q 
    
    res <- c(dS, dE, dI, dV, dQ)
    list(res)
  }
  )
}

SEIVQ_lograbies <- function(t, logx, params) { 
  ## Input:
  ##  - t - time
  ##  - logx - state variables
  ## Output: list of population sizes for each compartment of the model
  ## note: FD framework for the entire population
  
  # initialize
  x <- exp(logx)
  S <- x[1]
  E <- x[2]
  I <- x[3]
  V <- x[4]
  Q <- x[5]
  
  N <- S + E + I + V + Q
  
  # calculate population sizes
  with(as.list(params), { 
    dS <- b*N - d*S - delta*N*S - beta*S*I/N - vc*S + wn*V 
    dE <- beta*S*I/N - d*E - delta*N*E - sigma*E 
    dI <- sigma*E - d*I - delta*N*I - gamma*I - q*I
    dV <- vc*S - d*V - delta*N*V - wn*V 
    dQ <- q*I - d*Q - delta*N*Q - tau*Q 
    
    res <- c(dS/S, dE/E, dI/I, dV/V, dQ/Q)
    list(res)
  }
  )
}

## ## Code I (B) WORK ON THIS!
SEIRVQreactive_rabies <- function(t, x, params) { 
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
  
  # calculate population sizes
  with(as.list(params), {
    N <- S + E + I + V + Q
    q <- D
    
    dS <- b*N - d*S - delta*N*S - beta*S*I/N - vc*S + wn*V 
    dE <- beta*S*I/N - d*E - delta*N*E - sigma*E 
    dI <- sigma*E - d*I - delta*N*I - gamma*I - q*I
    dR <- gamma*I + tau*Q + d*S + d*E + d*I + d*V + d*Q + delta*N*S + delta*N*E + delta*N*I + delta*N*V + delta*N*Q
    dV <- vc*S - d*V - delta*N*V - wn*V 
    dQ <- q*I - d*Q - delta*N*Q - tau*Q 
    
    dD <- (omega*I-D)/eta
    
    res <- c(dS, dE, dI, dR, dV, dQ)
    list(res)
  }
  )
}

## Code II
# Newton's method to find equilibria of vector field.
# func() must have the same input arguments and returns as for lsoda/rk4.  
# x0 = intial guess at equilibrium. If x0 is not supplied in the call, 
# the user chooses it from the current graphics device via locator()
# and the equilibrium is plotted to the same device. Plotting
# symbol is closed/open=stable/unstable, circle/triangle=eigenvalues imaginary/real.   
# tol= Convergence tolerance 
# niter = Maximum number of iterations
# inc = finite-difference increment for derivative estimates 
# Coded 5/25/06 by SPE based on Matlab toggle.m by JG 
newton=function(func,x0=NULL,parms=NULL,tol=1e-14,niter=40,inc=1e-7) {
  x=x0; #initial x  
  if (is.null(x0)){x = locator(n=1); x=c(x$x,x$y)};
  nx = length(x); # length of state vector
  ######### Newton iteration loop: start  
  for(i in 1:niter){  
    y = func(0,x,parms)[[1]] 
    df = matrix(0,nx,nx); # Compute df
    for(j in 1:nx) {
      #Increment vector for estimating derivative wrt jth coordinate
      v=rep(0,nx); 
      v[j] = inc; 
      df[,j] = (func(t,x+v,parms)[[1]] - y)/inc;
    }
    if (sum(y^2) < tol){  #check for convergence 
      if(is.null(x0)){
        ev=eigen(df)$values; pch1=1+as.numeric(Im(ev[1])!=0); pch2=1+as.numeric(max(Re(ev))<0);
        pchs=matrix( c(2,17,1,16),2,2,byrow=T); 	
        points(x[1],x[2],type="p",pch=pchs[pch1,pch2],cex=1.5)
      }
      return(list(x=x,df=df))   
    }	
    x = x - solve(df,y) # one more step if needed 
    cat(i, x, "\n") #print out the next iterate 
  }
  ######### Newton iteration loop: end  
  cat("Convergence failed"); 
}

