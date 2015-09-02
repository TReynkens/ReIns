

# Hill estimator suitable for interval censored observations
# Z is the amount already paid and I are the incurred values
#
# if f(lower) and f(upper) have the same sign, we increase f(upper) at most maxup
# times by up

ciHill <- function(Z, I, censored, threshold = NULL, lower = 0.01, upper = 2,
                   up = 0.2, maxup = 50, warnings = FALSE, plot = FALSE, 
                   add = FALSE, main = "Estimates for the EVI", ...) {
  
  # Check input arguments
  checkInput(Z)
  
  if(length(censored)!=1) {
    if(length(Z) != length(censored)) {
      stop("Z and censored should have the same length.")
    }
  } else {
    censored <- rep(0,length(Z))
  }
  
  if(length(Z)!=length(I)) {
    stop("Z and I should have the same length.")
  } 
  
  Zsort <- sort(Z)
  n <- length(Z)

  # delta is 1 if a data point is not censored
  delta <- !(censored)
  
  null_thresh <- is.null(threshold)
  
  # Use order statistics when threshold=NULL
  if (null_thresh ) {
    K <- 1:(n-1)
    threshold <- Zsort[n-K]
  }
  
  gamma1 <- numeric(length(threshold))
  
  for (i in 1:length(threshold)) {
    ind <- which(Z>threshold[i])
    
    # Number of exceedances
    nt <- length(ind)
    
    # Proportion of non-censored over the treshold
    p <- 1/nt * sum(delta[ind])
    
    # Hill estimator adapted for censoring
    chill <- sum(log(Z[ind]/threshold[i]))/nt / p
    
    IZ <- I[ind]/Z[ind]

    # Function to compute roots of
    f__ci <- function(gamma1) {
      # Numerical problem when IZ=0, remove when not censored
      a <- ifelse(IZ==1 & !censored[ind], 1, 1-IZ^(-1/gamma1) )
      
      return( gamma1 + 1/p * 1/nt * sum( (1-delta[ind]) * (IZ)^(-1/gamma1) * log(IZ) / a) - chill )
    }

    gamma1[i] <- NA
    
    if (!is.na(f__ci(lower)) & !is.na(f__ci(upper)) ) {
      
      j <- 1
      # Increase upper if f(upper) and f(lower) have the same sign
      # No solution if f(upper) and f(lower) still have the same sign
      # in the end
      while ( (f__ci(upper) * f__ci(lower)) >= 0 & j<maxup) {
        j <- j+1
        upper <- upper + up
      }
      
      if ( (f__ci(upper) * f__ci(lower)) < 0) {
        gamma1[i] <- uniroot(f__ci, lower=lower, upper=upper)$root
      } else if(warnings) {
        warning("No upper found such that f(upper) and f(lower) have different signs.")
      }
      
    }

  }
  

  if (null_thresh) {
    # No threshold provided so order statistics are used.
    # Plots and output similar to other estimators

    plotfun(K, gamma1[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
    
    output(list(k=K, threshold=threshold, gamma1=gamma1), plot=plot, add=add)
    
  } else {
    # Threshold provided, only plots if more than 1 threshold is given
    
    if(length(threshold)>1 & (plot|add)) {
      # plots if TRUE
      plotfun(threshold, gamma1, type="l", xlab="Threshold", ylab="gamma1", main=main, plot=plot, add=add, ...)
    } else {
      plot <- add <- FALSE
    }
    output(list(threshold=threshold, gamma1=gamma1), plot=plot, add=add)
  }
  
}


# Weissman estimator for small exceedance probabilities for interval censored data
ciProb <- function(Z, I, censored, gamma1, q, plot = FALSE, add = FALSE, 
                      main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  checkInput(Z)
  censored <- checkCensored(censored, length(Z))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(Z) 
  n <- length(Zsort)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)
  
  # Estimator for probabilities
  prob[K] <- (1-tu) * (q/Zsort[n-K])^(-1/gamma1[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
}


# Weissman estimator for return periods for interval censored data
ciReturn <- function(Z, I, censored, gamma1, q, plot = FALSE, add = FALSE, 
                   main = "Estimates of return period", ...) {
  
  # Check input arguments
  checkInput(Z)
  censored <- checkCensored(censored, length(Z))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(Z) 
  n <- length(Zsort)
  R <- numeric(n)
  K <- 1:(n-1)
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)
  
  # Estimator for probabilities
  R[K] <- 1 / ( (1-tu) * (q/Zsort[n-K])^(-1/gamma1[K]) )
  R[R<1] <- NA
  
  # plots if TRUE
  plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  output(list(k=K, R=R[K], q=q), plot=plot, add=add)
}


#######################################################################################################

# Fit interval censored EPD to data using MLE.
# start is the starting value for optimisation: (gamma_start,kappa_start)
# delta is the non-censoring indicator
# notdelta is the censoring indicator, when NULL we set it to !delta

.ciEPDfit <- function(Zt, It, delta, notdelta = NULL, tau, start = c(0.1,1), warnings = FALSE) {
  
  if (is.numeric(start) & length(start)==2) {
    gamma_start <- start[1]
    kappa_start <- start[2]
  } else {
    stop("start should be a 2-dimensional numeric vector.")
  }
  
  delta <- as.logical(delta)
  if (is.null(notdelta)) {
    notdelta <- as.logical(!delta)
  }
  notdelta <- as.logical(notdelta)
  
  if(ifelse(length(Zt)>1,var(Zt)==0,0)) {
    sg <- c(NA,NA)
  } else {
    #Note that optim minimises a function so we use minus the log-likelihood function
    fit = optim(par=c(gamma_start,kappa_start), fn=.ciEPDneglogL, Zt=Zt, It=It, delta=delta, notdelta=notdelta, tau=tau)
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



# Minus log-likelihood for iid interval censored EPD random variables
# delta is the non-censoring indicator, make sure that delta
# and notdelta are logical!
.ciEPDneglogL <- function(theta, tau, Zt, It, delta, notdelta) {
  
  gamma <- theta[1]
  kappa <- theta[2]
  
  if (kappa<=max(-1,1/tau) | gamma<=0) {
    logL <- -10^6
  } else {

#     logL <- sum(delta * log(.dEPD(Zt, gamma=gamma, kappa=kappa, tau=tau)) ) +
#       sum( (1-delta) *  log(.pEPD(It, gamma=gamma, kappa=kappa, tau=tau) -  
#                   .pEPD(Zt, gamma=gamma, kappa=kappa, tau=tau)) )

    # Use delta and notdelta as index for the vector to avoid numerical issues
    logL <- sum(  log(.dEPD(Zt[delta], gamma=gamma, kappa=kappa, tau=tau)) ) +
            sum(  log(.pEPD(It[notdelta], gamma=gamma, kappa=kappa, tau=tau) -  
                             .pEPD(Zt[notdelta], gamma=gamma, kappa=kappa, tau=tau)) )
  }
  
  # minus log-likelihood for optimisation
  return(-logL)
}

# Interval censored EPD estimator for gamma_1 and kappa_1
#
# Z is the amount already paid and I are the incurred values
#
# rho is a parameter for the rho_estimator of Fraga Alves et al. (2003)
# when strictly positive or (a) choice(s) for rho if negative
#
# Start should be a 2-dimensional vector: (gamma_start,kappa_start)
# or NULL: we then use Hill estimates as initial values for gamma and kappa if threshold=NULL
# and c(1,0) otherwise

ciEPD <- function(Z, I, censored, threshold = NULL, rho = -1, beta = NULL, start = NULL, warnings = FALSE,
                  plot=FALSE, add=FALSE, main = "EPD estimates of EVI", ...){
                    
    
    # Check input arguments
    checkInput(Z)
    
    if(length(censored)!=1) {
      if(length(Z) != length(censored)) {
        stop("Z and censored should have the same length.")
      }
    } else {
      censored <- rep(0,length(Z))
    }
    
    if(length(Z)!=length(I)) {
      stop("Z and I should have the same length.")
    } 
    
    Zsort <- sort(Z)
    n <- length(Z)
    
    # delta is 1 if a data point is not censored
    delta <- !(censored)
    delta <- as.logical(delta)
    censored <- as.logical(censored)

    null_thresh <- is.null(threshold)
    
    HillZ <- Hill(Z)$gamma
    
    # Use order statistics when threshold=NULL
    if (null_thresh) {
      K <- 1:(n-1)
      threshold <- Zsort[n-K]
    }
    
    nrho <- length(rho)
    
    gamma1 <- matrix(0,length(threshold),nrho)
    kappa1 <- matrix(0,length(threshold),nrho)
    
    if (is.null(beta)) {
      # determine beta using rho
      
      beta <- matrix(0,n-1,nrho)
      
      if (all(rho>0) & nrho==1) {
        rho <- .rhoEst(Z,alpha=1,tau=rho)$rho
        
        # Estimates for rho of Fraga Alves et al. (2003) used 
        # and hence a different value of beta for each k
        beta[,1] <- -rho/HillZ
        
        for(j in 1:length(rho)) {
          beta[,j] <- -rho[j]/HillZ
        }
        
      } else if(all(rho<0)) {
        
        # rho is provided => beta is constant over k
        # (but differs with rho)
        for(j in 1:nrho) {
          beta[,j] <- -rho[j]/HillZ
        }
        
        
      } else {
        stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
      }
      
    } else {
      
      # use provided value(s) for beta
      if (length(beta)==1) {
        beta <- matrix(beta,n-1,nrho)
      } else if(length(beta)==n-1){
        beta <- as.matrix(rep(beta,nrho),ncol=length(rho),byrow=FALSE)
      } else {
        stop(paste0("beta should have length 1 or n-1 = ",n-1,"."))
      }
    }
    
    for(j in 1:length(rho)) {
      
      for(i in 1:length(threshold)) {
        
        # Starting values
        if(is.null(start)) {
          start2 <- numeric(2)
          if(null_thresh) {
            start2[1] <- HillZ[K[i]]
            start2[2] <- 0
          } else {
            start2[1] <- 1
            start2[2] <- 0
          }

        } else {
          start2 <- start
        }
        
        ind <- which(Z>threshold[i])
        
        Zt <- Z[ind]/threshold[i]
        It <- I[ind]/threshold[i]
        
        # tau is -beta
        res <- .ciEPDfit(Zt=Zt, It=It, tau=-beta[j], start=start2, delta=delta[ind],
                        notdelta=censored[ind], warnings=warnings)
        
        gamma1[i,j] <- res[1]
        kappa1[i,j] <- res[2]
      }
      
    }
    # Only positive values are allowed
    gamma1[gamma1<=0] <- NA
    
    

    
    if (null_thresh) {
      # No threshold provided so order statistics are used.
      # Plots and output similar to other estimators
      
      plotfun(K, gamma1[K,1], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
      if(length(rho)>1 & (plot | add)) {
        # Add lines
        for(j in 2:length(rho)) {
          lines(K,gamma1[K,j],lty=j)
        }
      }
      
    } else {
      # Threshold provided, only plots if more than 1 threshold is given
      
      if(length(threshold)>1 & (plot|add)) {
        # plots if TRUE
        plotfun(threshold, gamma1[,1], type="l", xlab="Threshold", ylab="gamma1", main=main, plot=plot, add=add, ...)
        if(length(rho)>1 & (plot | add)) {
          # Add lines
          for(j in 2:length(rho)) {
            lines(K,gamma1[,j],lty=j)
          }
        }
      } else {
        plot <- add <- FALSE
      }
    }
    
    
    # Transform to vectors if rho is a single value
    if(length(rho)==1) {
      gamma1 <- as.vector(gamma1)
      kappa1 <- as.vector(kappa1)
      beta <- as.vector(beta)
    }
    
    
    L <- list(gamma1=gamma1, kappa1=kappa1, beta=beta)
    if (null_thresh) {
      L$k <- K
    }
   
    output(L, plot=plot, add=add)
}


# Estimator for small exceedance probabilities for interval censored data
# using EPD estimates for interval censored data
ciProbEPD <- function(Z, I, censored, gamma1, kappa1, beta, q, plot = FALSE, add=FALSE,
                      main = "Estimates of small exceedance probability",...) {
             
  # Check input arguments
  checkInput(Z)
  censored <- checkCensored(censored, length(Z))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(Z) 
  n <- length(Zsort)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)                                  

  K2 <- K[!is.na(gamma1[K])]
  
  prob[K2] <- (1-tu) * (1-.pEPD(q/Zsort[n-K2],gamma=gamma1[K2],kappa=kappa1[K2],tau=-beta))
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}

# Estimator for sreturn periods for interval censored data
# using EPD estimates for interval censored data
ciReturnEPD <- function(Z, I, censored, gamma1, kappa1, beta, q, plot = FALSE, add=FALSE,
                      main = "Estimates of return period",...) {
  
  # Check input arguments
  checkInput(Z)
  censored <- checkCensored(censored, length(Z))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(Z) 
  n <- length(Zsort)
  R <- numeric(n)
  K <- 1:(n-1)
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)                                  
  
  K2 <- K[!is.na(gamma1[K])]
  
  R[K2] <- 1 / ((1-tu) * (1-.pEPD(q/Zsort[n-K2],gamma=gamma1[K2],kappa=kappa1[K2],tau=-beta)) )
  R[R<1] <- NA
  
  # plots if TRUE
  plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  output(list(k=K, R=R[K], q=q), plot=plot, add=add)
  
}

