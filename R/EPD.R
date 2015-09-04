
########################################################
# Distribution

# Check input for EPD distribution
.EPDinput <- function(y, gamma, kappa, tau, kappaTau = TRUE) {
  
  # Check if arguments are numeric
  if(!is.numeric(gamma)) {
    stop("gamma should be numeric.")
  }
  
  if (!is.numeric(kappa)) {
    stop("kappa should be numeric.")
  }
  
  if (!is.numeric(tau)) {
    stop("tau should be numeric.")
  }
  
  # Check if right sign
  if (any(tau>=0)) {
    stop("tau should be strictly negative.")
  }
  
  if (any(gamma<=0)) {
    stop("gamma should be strictly positive.")
  }
  
  if (kappaTau) {
    if (any(kappa<=pmax(-1,1/tau))) {
      stop("kappa should be larger than max(-1,1/tau).")
    }
  }
  
  # Check if correct length
  ly <- length(y)
  lg <- length(gamma)
  lk <- length(kappa)
  lt <- length(tau)
  
  l <- c(ly, lg, lk, lt)
  
  # Indices of lengths larger than 1
  ind <- which(l>1)
  
  if (length(ind)>1) {
    # Check that lengths larger than 1 are equal
    if(!all.equal(l[ind])) {
      stop("All input arguments should have length 1 or equal length.")
    }
  }
  

}
# Density of an extended Pareto distribution
.dEPD <- function(x, gamma, kappa, tau = -1, log = FALSE) {
  
  # Check input
  .EPDinput(x, gamma, kappa, tau, kappaTau = TRUE)
  
  # Compute density
  d <- 1 / (gamma*x^(1/gamma+1)) * (1+kappa*(1-x^tau))^(-1/gamma-1) * 
    (1+kappa*(1-(1+tau)*x^tau))
  # Formula is not valid for values below 1
  d[x<=1] <- 0
  
  if (log) d <- log(d)
  
  return(d)
}

# CDF of an extended Pareto distribution
.pEPD <- function(x, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
  
  # Check input
  .EPDinput(x, gamma, kappa, tau, kappaTau = FALSE)

  # Compute probabilities
  p <- 1 - (x * (1+kappa*(1-x^tau)))^(-1/gamma) 
  # Formula is not valid for values below 1
  p[x<=1] <- 0
  
  # Problems when condition not satisfied
  if (any(kappa<=pmax(-1,1/tau))) {
    if (length(kappa)>1 | length(tau)>1) {
      p[kappa<=pmax(-1,1/tau)] <- NA
    } else {
      p <- NA
    }
  }
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}

# Quantile function of EPD
.qEPD <-  function(p, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
  
  # Check input
  .EPDinput(p, gamma, kappa, tau, kappaTau = TRUE)
  

  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  if (any(p<0 | p>1)) {
    stop("p should be between 0 and 1.")
  }
  
  # Compute quantiles numerically
  l <- length(p)
  Q <- numeric(l)
  
  # Take 100 as endpoint for interval to search over unless not large enough
  endpoint <- ifelse(.pEPD(100,gamma,kappa,tau)>=max(p[p<1]), 100, 1000)
  
  for(i in 1:l) {
    
    
    if (p[i]<10^(-14)) {
      # p=0 case
      Q[i] <- 1
      
    } else if (p[i]<1) {
      # 0<p<1 case
      
      # Function to compute roots of
      f <- function(x) {
        (1-p[i])^(-gamma) - x*(1+kappa*(1-x^tau))
      }
      # If root solving fails return NA
      Q[i] <- tryCatch(uniroot(f,lower=1,upper=endpoint)$root, error=function(e) NA) 
      
    } else {
      # p=1 case
      Q[i] <- Inf
    }
    
  }
  
  return(Q)
}



# Random number generation for EPD
.rEPD <-  function(n, gamma, kappa, tau = -1) {

  # Rely on input checking in .qEPD
  
  # Generate random numbers
  return(.qEPD(runif(n), gamma=gamma, kappa=kappa, tau=tau))
}

########################################################
# Estimation


# Start should be a 2-dimensional vector: (gamma_start,kappa_start)
# or NULL (we then use Hill estimates as initial values for gamma and kappa)
# It is only used when direct=TRUE.
#
# rho is a parameter for the rho_estimator of Fraga Alves et al. (2003)
# when strictly positive or (a) choice(s) for rho if negative
EPD <- function(data, rho = -1, start = NULL, direct = FALSE, warnings = FALSE, 
                plot = FALSE, add = FALSE, main = "EPD estimates of EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  if (direct) {
    # Select parameters by minimising MLE
    EPD <- .EPDdirectMLE(data=data, rho=rho, start=start, warnings=warnings)
  } else {
    # Select parameter using approach of Beirlant, Joosens and Segers (2009). 
    EPD <- .EPDredMLE(data=data, rho=rho)
  }
  
  # plots if TRUE
  .plotfun(K, EPD$gamma[K,1], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  
  # Transform to vectors if rho is a single value
  if(length(rho)==1) {
    EPD$gamma <- as.vector(EPD$gamma)
    EPD$kappa <- as.vector(EPD$kappa)
    EPD$tau <- as.vector(EPD$tau)
    
  } else if(plot | add) {
    # Add lines
    for(j in 2:length(rho)) {
      lines(K,EPD$gamma[K,j],lty=j)
    }
  }
  
  if(length(rho)==1) {
    .output(list(k=K, gamma=EPD$gamma[K], kappa=EPD$kappa[K], tau=EPD$tau[K]), plot=plot, add=add)
  } else {
    .output(list(k=K, gamma=EPD$gamma[K,], kappa=EPD$kappa[K,], tau=EPD$tau[K,]), plot=plot, add=add)
  }
}


# Fit EPD using approach of Beirlant, Joosens and Segers (2009). 
.EPDredMLE <- function(data, rho = -1) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  nrho <- length(rho)
  rho.orig <- rho
  
  H <- Hill(data, plot=FALSE)$gamma
  
  if (all(rho>0) & nrho==1) {
    rho <- .rhoEst(data,alpha=1,tau=rho)$rho
    beta <- -rho

  } else if(all(rho<0)) {
    beta <- -rho
    
  } else {
    stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
  }
  
  
  gamma <- matrix(0,n-1,nrho)
  kappa <- matrix(0,n-1,nrho)
  tau <- matrix(0,n-1,nrho)
  
  beta <- -rho
  
  for(j in 1:nrho) {
    
    # tau
    if (nrho==1 & all(rho.orig>0)) {
      # Estimates for rho of Fraga Alves et al. (2003) used 
      # and hence a different value of beta for each k
      tau[K,1] <- -beta[K]/H[K]
      
      # rho differs with k
      rhovec <- rho
    } else {
      # rho is provided => beta is constant over k
      # (but differs with rho)
      tau[K,j] <- -beta[j]/H[K]
      
      # rho is constant over k
      rhovec <- rep(rho[j],n-1)
    }

    # kappa
    E <- numeric(n-1)
    
    for(k in K) {
      i <- 1:k
      E[k] <- 1/k * sum( (X[n-k+i]/X[n-k])^tau[k,j] )
    }
    
    kappa[K,j] <- H[K] * (1-2*rhovec[K]) * (1-rhovec[K])^3 / rhovec[K]^4 * (E - 1 / (1-rhovec[K]))

    
    # gamma
    gamma[K,j] <- H[K] - kappa[K,j] * rhovec[K] / (1 - rhovec[K])
      
  }
  
  return(list(gamma=gamma, kappa=kappa, tau=tau))
}




# Fit EPD by minimising MLE
.EPDdirectMLE <- function(data, rho = -1, start = NULL,  warnings = FALSE) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)

  if (n==1) {
    stop("We need at least two data points.")
  }
  
  nrho <- length(rho)
  rho.orig <- rho
  
  H <- Hill(data, plot=FALSE)$gamma
  
  if (all(rho>0) & nrho==1) {
    rho <- .rhoEst(data,alpha=1,tau=rho)$rho
    beta <- -rho
    
  } else if(all(rho<0)) {
    beta <- -rho
    
  } else {
    stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
  }
  
  gamma <- matrix(0,n-1,nrho)
  kappa <- matrix(0,n-1,nrho)
  tau <- matrix(0,n-1,nrho)
  
  for(j in 1:nrho) {

    # Compute gamma and kappa for several values of k
    for(k in (n-1):1) {
      epddata <- data[data>X[n-k]]/X[n-k]
      
      # tau
      if (nrho==1 & all(rho.orig>0)) {
        # Estimates for rho of Fraga Alves et al. (2003) used 
        # and hence a different value of beta for each k
        tau[k,1] <- -beta[k]/H[k]
        
      } else {
        # rho is provided => beta is constant over k
        # (but differs with rho)
        tau[k,j] <- -beta[j]/H[k]
      }
      
      if(is.null(start)) {
        start2 <- numeric(2)
        start2[1] <- H[k]
        start2[2] <- 0
      } else if(is.matrix(start)) {
        
        if(nrow(start>=n-1)) {
          start2 <- numeric(2)
          start2[1] <- start[k,1]
          start2[2] <- start[k,2]
        } else {
          stop("start does not contain enough rows.")
        }
        
      } else {
        start2 <- start
      }
      
      if (tau[k,j]<0) {
        tmp <- EPDfit(epddata, start=start2, tau=tau[k,j])
        gamma[k,j] <- tmp[1]
        kappa[k,j] <- tmp[2]
      } else {
        # Problems if tau is >0
        gamma[k,j] <- kappa[k,j] <- NA
      }

    }
    
  }
  
  return(list(gamma=gamma,kappa=kappa))
}



# Fit EPD to data using MLE.
# start is the starting value for optimisation: (gamma_start,kappa_start)
EPDfit <- function(data, tau, start = c(0.1,1), warnings = FALSE) {
  
  if (is.numeric(start) & length(start)==2) {
    gamma_start <- start[1]
    kappa_start <- start[2]
  } else {
    stop("start should be a 2-dimensional numeric vector.")
  }
  
  
  if(ifelse(length(data)>1,var(data)==0,0)) {
    sg <- c(NA,NA)
  } else {
    
    #Note that optim minimises a function so we use minus the log-likelihood function
    fit = optim(par=c(gamma_start,kappa_start), fn=.EPDneglogL, Y=data, tau=tau)
    # fit = nlminb(start=c(gamma_start,log(sigma_start)),objective=neglogL, Y=data)
    sg = fit$par
    
    if (fit$convergence>0 & warnings) {
      warning("Optimisation did not complete succesfully.")
      if (!is.null(fit$message)) {
        print(fit$message)
      }
    }
  }
  return(sg)
}



# Minus log-likelihood for a (univariate sample Y) of iid EPD random variables
.EPDneglogL <- function(theta, Y, tau) {
  
  gamma <- theta[1]
  # Makes sure that sigma is positive
  kappa <- theta[2]

  if (kappa<=max(-1,1/tau) | gamma<=0) {
    logL <- -10^6
  } else {
    logL <- sum( log(.dEPD(Y, gamma=gamma, kappa=kappa, tau=tau)) )
  }
  
  # minus log-likelihood for optimisation
  return(-logL)
}

# Estimator for rho of Fraga Alves, Gomes and de Haan (2003)
.rhoEst <- function(data, alpha = 1, theta1 = 2, theta2 = 3, tau = 1) {
  
  # Check input arguments
  .checkInput(data)
  
  if(alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  
  if(tau<=0) {
    stop("tau should be strictly positive.")
  }
  
   
  X <- as.numeric(sort(data))
  n <- length(X)
  rho <- numeric(n)
  Tn <- numeric(n)
  K <- 1:(n-1)
  
  M_alpha <- numeric(n)
  M_alpha_theta1 <- numeric(n)
  M_alpha_theta2 <- numeric(n)
  
  l <- log(X[n-K+1])
  for(k in K) {
    M_alpha[k] <- sum( (l[1:k]-log(X[n-k]))^alpha ) / k
    M_alpha_theta1[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta1) ) / k
    M_alpha_theta2[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta2) ) / k
  }
  
  Tn[K] <- ( (M_alpha[K]/gamma(alpha+1))^tau - (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1)  ) / 
           ( (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1) - (M_alpha_theta2[K]/gamma(alpha*theta2+1))^(tau/theta2)  ) 
  
  rho[K] <- 1 - ( 2 * Tn[K] / ( 3 - Tn[K]) ) ^ (1/alpha)
   
  return(list(k=K, rho=rho[K], Tn=Tn[K]))
}

############################################################
# Probabilities and return periods


ProbEPD <- function(data, q, gamma, kappa, tau, plot = FALSE, add = FALSE,
                    main = "Estimates of small exceedance probability",...) {
  # Check input arguments
  .checkInput(data)
  
  if( length(gamma)!=length(kappa) | length(gamma)!=length(tau)) {
    stop("gamma, kappa and tau should have equal length.")
  } 
  
  X <- as.numeric(sort(data))
  n <- length(X)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  K2 <- K[which(gamma[K]>0)]
  
  prob[K2] <- (K2+1)/(n+1) * (1 - .pEPD(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}



ReturnEPD <- function(data, q, gamma, kappa, tau, plot = FALSE, add = FALSE, 
                      main = "Estimates of return period", ...) {
  # Check input arguments
  .checkInput(data)
  
  if( length(gamma)!=length(kappa) | length(gamma)!=length(tau)) {
    stop("gamma, kappa and tau should have equal length.")
  } 
  
  X <- as.numeric(sort(data))
  n <- length(X)
  r <- numeric(n)
  K <- 1:(n-1)
  
  K2 <- K[which(gamma[K]>0)]
  
  r[K2] <- (n+1)/(K2+1) / (1 - .pEPD(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))
  
  
  r[which(gamma[K]<=0)] <- NA
  
  r[r<=0] = NA
  
  
  # plots if TRUE
  .plotfun(K, r[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=r[K], q=q), plot=plot, add=add)
  
}


