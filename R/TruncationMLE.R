# Estimate gamma and tau using MLE for truncated model
# eps is the numerical tolerance used in the likelihood
trMLE <- function(data, start = c(1,1), eps = 10^(-10), 
                  plot = TRUE, add = FALSE, main = "Estimates for EVI", ...) {
  
  X <- sort(data)
  n <- length(X)
  gamma <- numeric(n)
  tau <- numeric(n)
  beta <- numeric(n)
  K <- 1:(n-1)
  
  # Convergence indicator
  conv <- numeric(n)
  
  start_orig <- start

  # Start in back to use previous estimates as starting values
  for (k in (n-1):2) {
    
    # Use provided starting value if k=n-1
    if (k!=(n-1)) {
      start <- c(gamma[k+1],tau[k+1])
    }
    
    E <- X[n-(1:k)+1] - X[n-k]
    
    E[1] <- X[n] - X[n-k]
    
    # Minus log-likelihood
    lik <- function(x) {
      gammapar <- x[1]
      taupar <- x[2]
      
      a <- 1+taupar*E[-1]
      beta <- 1-(1+taupar*E[1])^(-1/gammapar)
      
      # tau and gamma need to have the same sign
      # and a cannot be <=0
      if (gammapar/taupar<eps | min(a)<eps | 1+taupar*E[1]<eps) {
        L <- -10^6
        
      } else {
        
        L <- (k-1)*log(taupar/gammapar) - (1+1/gammapar)*sum(log(a)) - (k-1)*log(beta)
      }
      return(-L)
    }
    
    # Check if suitable starting value
    if (!is.numeric(lik(start)) | !is.finite(lik(start)) | is.nan(lik(start))
        | lik(start)==10^6) {
      start <- start_orig
    }
    
    # Minimise minus log-likelihood
    sol <- optim(start, fn=lik, method = "Nelder-Mead")
    
    # Check convergence
    conv[k] <- sol$conv
    
    # Obtained estimate for gamma
    gamma[k] <- sol$par[1]
    
    # Obtained estimate for tau
    tau[k] <- sol$par[2]
  }
  
  
  ### plots if TRUE  
  .plotfun(K, gamma[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  
  
  .output(list(k=K, gamma=gamma[K], tau=tau[K], sigma=gamma[K]/tau[K], conv=conv), plot=plot, add=add)
}



# Estimator for odds ratio
trDTMLE <- function(data, gamma, tau, plot = FALSE, add = FALSE, main = "Estimates of DT", ...) {
  
  X <- sort(data)
  n <- length(X)
  K <- 1:(n-1)
  
  E <- X[n] - X[n-K]
  
  # Formula 17
  DT <- pmax(0, K/n * ((1+tau[K]*E)^(-1/gamma[K])-1/K) / (1-(1+tau[K]*E)^(-1/gamma[K])))
  
  ### plots if TRUE  
  .plotfun(K, DT[K], type="l", xlab="k", ylab=expression(D[T]), 
           main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, DT=DT[K]), plot=plot, add=add)
  
  
}

# Estimates of quantile Q(1-p)
trQuantMLE <- function(data, gamma, tau, DT, p, Y = FALSE, plot = FALSE, add = FALSE, main = "Estimates of extreme quantile", ...) {
  
  X <- sort(data)
  n <- length(X)
  K <- 1:(n-1)
  
  if (Y) {
    
    # Quantile of Y, formula 19
    Q <- X[n-K] + 1/tau[K] * (((DT[K]+K/n) / (p*(DT[K]+1)))^(gamma[K]) - 1)
    
  } else {
    
    # Quantile of X, formula 18
    Q <- X[n-K] + 1/tau[K] * (((DT[K]+K/n) / (DT[K]+p))^(gamma[K]) - 1)
    
  }
  
  
  ### plots if TRUE  
  .plotfun(K, Q[K], type="l", xlab="k", ylab=ifelse(Y, expression(Q[Y](1-p)), expression(Q(1-p))), 
           main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, q=Q[K]), plot=plot, add=add)
}


# Estimates of endpoint
trEndpointMLE <- function(data, gamma, tau, plot = FALSE, add = FALSE, main = "Estimates of endpoint", ...) {
  
  X <- sort(data)
  n <- length(X)
  K <- 1:(n-1)
  
  
  a <- (1+tau[K]*(X[n]-X[n-K]))^(-1/gamma[K])
  # Quantile of X, formula 18
  Tn <- X[n-K] + 1/tau[K] * (((1-1/K) / (a-1/K))^(gamma[K]) - 1)

  ### plots if TRUE  
  .plotfun(K, Tn[K], type="l", xlab="k", ylab=expression(T[k]), 
           main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, Tk=Tn[K]), plot=plot, add=add)
}


### exceedance probability estimation
trProbMLE <- function(data, gamma, tau, DT, q, plot = FALSE, add = FALSE,
                      main="Estimates of small exceedance probability", ...) {
  
  
  X <- sort(data)
  n <- length(X)
  K <- 1:(n-1)
  prob <- numeric(n)
  
  
  # Formula 22
  prob[K] <- (1+DT[K]) * K/n * (1+tau[K]*(q-X[n-K]))^(-1/gamma[K]) - DT[K]
  
  
  
  ### plots if TRUE
  
  if( !all(is.na(prob[K])) ) {
    .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  }
  
  
  ### output list with values of k, exceedance probabilitye estimates 
  ### and the considered high quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}

# Test for truncation
trTestMLE <- function(data, gamma, tau, alpha = 0.05, plot = TRUE, main = "Test for truncation", ...) {
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 2:(n - 1)
  
  E1k <- X[n] - X[n-K]
  tv <- K * (1 + tau[K] * E1k)^(-1/gamma[K])
  cv <- log(1 / alpha)
  reject <- (tv > cv)
  Pval <- exp(-tv)
  #Pval <- (1-tv/K)^K
  if (plot) {
    plot(K, Pval[K], ylim = c(0, 1), type = "l", xlab = "k", 
         ylab = "P-value", main = main, ...)
    abline(h = alpha, col = "blue", ...)
  }
  
  return(list(k = K, testVal = tv, critVal = cv, Pval = Pval, reject = reject))
}
