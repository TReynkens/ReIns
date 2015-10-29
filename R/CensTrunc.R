
#Hill estimator suitable for censored and truncated observations
#Z is the amount already paid and I are the incurred values
trcHill <- function(data, censored, r = 1, threshold = NULL, lower = 0.01, upper = 2, 
                    up = 0.2, maxup = 50, plot = FALSE, add = FALSE, 
                     main = "Estimates for the EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  if(length(censored)!=1) {
    if(length(data) != length(censored)) {
      stop("data and censored should have the same length.")
    }
  } else {
    censored <- rep(0,length(data))
  }
  

  Zsort <- sort(data)
  n <- length(Zsort)
  
  # delta is 1 if a data point is not censored
  delta <- !(censored)
  
  null_thresh <- is.null(threshold)
  
  # Use order statistics when threshold=NULL
  if (null_thresh ) {
    K <- r:(n-1)
    threshold <- Zsort[n-K]
  }
  
  
  gamma1 <- numeric(length(threshold))
  
  for (i in 1:length(threshold)) {
    
    # Index for data which are larger than t
    ind <- which(data>threshold[i])
    
    # Number of exceedances
    nt <- length(ind)
    
    # Proportion of non-censored over the treshold
    p <- 1/nt * sum(delta[ind])
    
    # Hill estimator adapted for censoring
    chill <- sum(log(data[ind]/threshold[i]))/nt / p
    
    beta <- Zsort[n]/threshold[i]
    
  
#     # Function to compute roots of
#     f_trc <- function(gamma1) {
#       gamma1 + 1/p * (log(1/beta) * beta^(-1/gamma1)) / (1-beta^(-1/gamma1)) - chill
#     }

    Zt <- data[ind]/threshold[i]
    
    # Function to compute roots of
    f_trc <- function(gamma1) {
      # Numerical problem when Zt=beta, contribution removed then since log(1)=0
      a <- Zt^(-1/gamma1) - beta^(-1/gamma1)
      a[abs(Zt-beta)<10^(-8)] <- 1
      
      return( gamma1 - 1/(nt*p) * sum( (1-delta[ind]) * beta^(-1/gamma1) * log(Zt/beta) / a ) 
              + 1/p * (log(1/beta) * beta^(-1/gamma1)) / (1-beta^(-1/gamma1)) - chill )
    }

    
    gamma1[i] <- 10^(-8)
    
    if (!is.na(f_trc(lower)) & !is.na(f_trc(upper)) ) {
      
      j <- 1
      # Increase upper if f(upper) and f(lower) have the same sign
      # No solution if f(upper) and f(lower) still have the same sign
      # in the end
      while ( (f_trc(upper) * f_trc(lower)) >= 0 & j<maxup) {
        j <- j+1
        upper <- upper + up
      }
      
      if ( (f_trc(upper) * f_trc(lower)) < 0) {
        gamma1[i] <- uniroot(f_trc, lower=lower, upper=upper)$root
      } 
      
    }
    
  }
  
  
  if (null_thresh) {
    # No threshold provided so order statistics are used.
    # Plots and output similar to other estimators
    
    .plotfun(K, gamma1[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
    
    .output(list(k=K, threshold=threshold, gamma1=gamma1), plot=plot, add=add)
    
  } else {
    # Threshold provided, only plots if more than 1 threshold is given
    
    if(length(threshold)>1 & (plot|add)) {
      # plots if TRUE
      .plotfun(threshold, gamma1, type="l", xlab="Threshold", ylab="gamma1", main=main, plot=plot, add=add, ...)
    } else {
      plot <- add <- FALSE
    }
    .output(list(threshold=threshold, gamma1=gamma1), plot=plot, add=add)
  }
  
}





#Hill estimator suitable for interval censored and truncated observations
#Z is the amount already paid and I are the incurred values
trciHill <- function(Z, I, censored, r = 1, threshold = NULL, lower = 0.01, upper = 2, 
                      up = 0.2, maxup = 50, warnings = FALSE, plot = FALSE, add = FALSE, 
                      main = "Estimates for the EVI", ...) {
  
  # Check input arguments
  .checkInput(Z)
  
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
    K <- r:(n-1)
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
    
    beta <- Zsort[n]/threshold[i]
    
    # Function to compute roots of
    f_trci <- function(gamma1) {
      # Numerical problem when IZ=1, contribution removed then since log(1)=0
      a <- 1-IZ^(-1/gamma1)
      a[IZ==1] <- 1
      
      return( gamma1 + 1/p * 1/nt * sum( (1-delta[ind]) * (IZ)^(-1/gamma1) * log(IZ) / a ) +
        1/p * (log(1/beta) * beta^(-1/gamma1)) / (1-beta^(-1/gamma1)) - chill )
    }
    
    gamma1[i] <- 10^(-8)
    
    if (!is.na(f_trci(lower)) & !is.na(f_trci(upper)) ) {
      
      j <- 1
      # Increase upper if f(upper) and f(lower) have the same sign
      # No solution if f(upper) and f(lower) still have the same sign
      # in the end
      while ( (f_trci(upper) * f_trci(lower)) >= 0 & j<maxup) {
        j <- j+1
        upper <- upper + up
      }
      
      if ( (f_trci(upper) * f_trci(lower)) < 0) {
        gamma1[i] <- uniroot(f_trci, lower=lower, upper=upper)$root
      } 
      
    }
    
  }
  
  
  if (null_thresh) {
    # No threshold provided so order statistics are used.
    # Plots and output similar to other estimators
    
    .plotfun(K, gamma1[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
    
    .output(list(k=K, threshold=threshold, gamma1=gamma1), plot=plot, add=add)
    
  } else {
    # Threshold provided, only plots if more than 1 threshold is given
    
    if(length(threshold)>1 & (plot|add)) {
      # plots if TRUE
      .plotfun(threshold, gamma1, type="l", xlab="Threshold", ylab="gamma1", main=main, plot=plot, add=add, ...)
    } else {
      plot <- add <- FALSE
    }
    .output(list(threshold=threshold, gamma1=gamma1), plot=plot, add=add)
  }
  
}


# Weissman estimator for small exceedance probabilities for censored and truncated data
trcProb <- function(data, censored, gamma1, q, r = 1, plot = FALSE, add = FALSE, 
                        main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(data) 
  n <- length(Zsort)
  prob <- numeric(n)
  K <- r:(n-1)
  
  # Kaplan-Meier estimator for survival function in Zsort[n-K]
  km <- KaplanMeier(Zsort[n-K], data = data, censored = censored)
  
  beta <- Zsort[n] / Zsort[n-K]
   
  # Estimator for probabilities
  prob[K] <- (1-km) * ( (q/Zsort[n-K])^(-1/gamma1[K]) - beta^(-1/gamma1[K]) ) / (1 - beta^(-1/gamma1[K]))
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
}


# Weissman estimator for small exceedance probabilities for interval censored and truncated data
trciProb <- function(Z, I, censored, gamma1, q, r = 1, plot = FALSE, add = FALSE, 
                           main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(Z)
  censored <- .checkCensored(censored, length(Z))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  Zsort <- sort(Z) 
  n <- length(Zsort)
  prob <- numeric(n)
  K <- r:(n-1)

  if(q>=Zsort[n]) {
    warning("q cannot be larger than the maximum of the data.")
  }
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)$cdf
  
  beta <- Zsort[n] / Zsort[n-K]
  
  # Estimator for probabilities
  prob[K] <- (1-tu) * ( (q/Zsort[n-K])^(-1/gamma1[K]) - beta^(-1/gamma1[K]) ) / (1 - beta^(-1/gamma1[K]))
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
}


trciDT <- function(Z, I, censored, gamma1, r = 1, plot=FALSE, add=FALSE, main="Estimates of DT", ...) {
  
  # Check input arguments
  .checkInput(Z)
  censored <- .checkCensored(censored, length(Z))
  
  Zsort <- as.numeric(sort(Z))
  n <- length(Z)
  DT <- numeric(n)
  K <- r:(n-1)
  
  if(length(gamma1)!=length(K)) {
    stop(paste("gamma1 should have length", length(K)))
  }
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)$cdf
  
  R <- Zsort[n-K] / Zsort[n-r+1]
  
  A <- 1/gamma1
  R <- Zsort[n-K]/Zsort[n-r+1]
  
  DT[K] <-  pmax( (1-tu) * (R^A-r/(K+1)) / (1-R^A), 0)
  
  ### plots if TRUE  
  .plotfun(K, DT[K], type="l", xlab="k", ylab="DT", main=main, plot=plot, add=add, ...)
  
  
  ### output list with values of k and
  ### corresponding estimates for DT
  
  .output(list(k=K, DT=DT[K]),plot=plot,add=add)
  
}


trciEndpoint <- function(Z, I, censored, gamma1, DT, r = 1, plot = FALSE, add = FALSE, 
                       main = "Estimates of Endpoint", ...) {
  
  # Check input arguments
  .checkInput(Z,gamma=gamma1,DT=DT,r=r)
  censored <- .checkCensored(censored, length(Z))
  
  Zsort <- as.numeric(sort(Z))
  n <- length(Z)
  Tk <- numeric(n)
  K <- r:(n-1)
  
  
  
  if(length(gamma1)!=length(K)) {
    stop(paste("gamma1 should have length", length(K)))
  }
  
  if(length(DT)!=length(K)) {
    stop(paste("DT should have length", length(K)))
  }
  
  # Turnbull estimator for CDF in Zsort[n-K]
  tu <- Turnbull(Zsort[n-K], L=Z, R=I, censored=censored)$cdf
  
  
  
  Tk[K] <-  exp( pmax( log(Zsort[n-K]) + gamma1 * log(1+(1-tu)/DT), log(Zsort[n])) )
  
  ### plots if TRUE  
  .plotfun(K, Tk[K], type="l", xlab="k", ylab="Tk", main=main, plot=plot, add=add, ...)
  
  ### output list with values of k and
  ### corresponding estimates for DT
  
  .output(list(k=K, Tk=Tk[K]),plot=plot,add=add)
  
}