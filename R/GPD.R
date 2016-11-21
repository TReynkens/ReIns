
# Computes the Peaks-Over-Threshold estimates of gamma using ML (Section 5.3.2 in Beirlant et al. (2004)) 
# Based on POT_KH.R (written by Klaus Herrmann) and POT.R.

###########################################################################

# minus log-likelihood for a (univariate sample Y) of iid GP random variables
.POTneglogL <- function(theta, Y) {
  
  gamma <- theta[1]
  # Makes sure that sigma is positive
  sigma <- exp(theta[2])
  
  n <- length(Y)
  yvec <- Y
    
  if (abs(gamma) < .Machine$double.eps) {
    logL <- -n*theta[2]-(1/sigma)*sum(Y)
  } else {
    yvec <- 1+gamma/sigma*Y
    
    if (min(yvec) > 0) {
      logL <- -n*theta[2]-(1/gamma+1)*sum(log(yvec))
    } else {
      logL <- -10^10   #Very small value  
    } 
  }

  # minus log-likelihood for optimisation
  return(-logL)
}

# Fit GPD to data using MLE.
# start is the starting value for optimisation: (gamma_start,sigma_start)
GPDfit <- function(data, start = c(0.1, 1), warnings = FALSE) {
  
  if (is.numeric(start) & length(start)==2) {
    if (start[2] <= 0) {
      stop("The starting value for sigma should be >0.")
    }
    gamma_start <- start[1]
    sigma_start <- start[2]
  } else {
    stop("start should be a 2-dimensional numeric vector.")
  }
  

  if (ifelse(length(data)>1, var(data)==0, 0)) {
    sg <- c(NA, NA)
  } else {
    
    #Note that optim minimises a function so we use minus the log-likelihood function
    fit <- optim(par=c(gamma_start, log(sigma_start)), fn=.POTneglogL, Y=data)
    # fit = nlminb(start=c(gamma_start,log(sigma_start)),objective=neglogL, Y=data)
    sg <- fit$par
    if (fit$convergence>0 & warnings) {
      warning("Optimisation did not complete succesfully.")
      if (!is.null(fit$message)) {
        print(fit$message)
      }
   }

    #We use log(sigma) in neglogL
    sg[2] <- exp(sg[2]) 
     
  }
  return(sg)
}


# Start should be a 2-dimensional vector: (gamma_start,sigma_start)
GPDmle <- function(data, start = c(0.1, 1), warnings = FALSE, logk = FALSE, 
                   plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  POT <- matrix(0, n, 2)
  K <- 1:(n-1)
  
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  #Compute gamma and sigma for several values of k
  for(k in (n-1):1) {
    potdata <- data[data>X[n-k]]-X[n-k]
    if (length(potdata)==0) {
      POT[k,] <- NA
    } else {
      POT[k,] <- GPDfit(potdata, start=start)
    }

  }

  # plots if TRUE
  if (logk) {
    .plotfun(log(K), POT[K,1], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, POT[K,1], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }
  
  .output(list(k=K, gamma=POT[K,1], sigma=POT[K,2]), plot=plot, add=add)
}


POT <- GPDmle






#####################################################################################################



# Computes estimates of small exceedance probability 1-F(q) 
# using POT estimators for gamma and sigma

ProbGPD <- function(data, gamma, sigma, q, plot = FALSE, add = FALSE, 
                    main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data, gamma, gammapos=FALSE)
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Estimator for tail probabilities

  # Formula 5.26
  prob[K] <- ((K+1)/(n+1)) * (1+gamma[K]/sigma[K]*(q-X[n-K]))^(-1/gamma[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}


ReturnGPD <- function(data, gamma, sigma, q, plot = FALSE, add = FALSE, 
                    main = "Estimates of large return period", ...) {
                      
    # Check input arguments
    .checkInput(data, gamma, gammapos=FALSE)
    
    if (length(q)>1) {
      stop("q should be a numeric of length 1.")
    }
    
    X <- as.numeric(sort(data))
    n <- length(X)
    r <- numeric(n)
    K <- 1:(n-1)
    
    # Estimator for tail probabilities
    
    # Formula 5.26
    r[K] <- (n+1)/(K+1) * (1+gamma[K]/sigma[K]*(q-X[n-K]))^(1/gamma[K])
    r[r<=0] <- NA
    
    # plots if TRUE
    .plotfun(K, r[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
    
    # output list with values of k, corresponding return period estimates 
    # and the considered large quantile q
    
    .output(list(k=K, R=r[K], q=q), plot=plot, add=add)
    
  }


##############################################################################
# Computes estimates of extreme quantile Q(1-p) 
# (Section 5.5.1 in Beirlant et al. (2004)) for a numeric vector of observations 
# (data) and as a function of k
#
# Estimates are based on prior estimates gamma for the
# extreme value index (Hill)
#
# If aest=1, estimates of a are based on (2.15) - gamma genHill
# or genZipf, if aest=2, estimates are based on (5.24) - gamma ERM
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot


QuantGPD <- function(data, gamma, sigma, p, plot = FALSE, add = FALSE, 
                     main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data, gamma, gammapos=FALSE)
  
  .checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  quant <- numeric(n)
  K <- 1:(n-1)
    
  # Estimator for extreme quantiles

  # Formula 5.23
  quant[K] <- X[n-K] + sigma[K] / gamma[K] * ( ((K+1)/((n+1)*p))^gamma[K] -1)
  
  # plots if TRUE
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p), plot=plot, add=add)
}

##############################################################################

# GPD residual plot when applying POT with threshold t
GPDresiduals <- function(data, t, gamma, sigma, plot = TRUE, main = "GPD residual plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  
  if (length(gamma)!=1 | length(sigma)!=1) {
    stop("gamma and sigma should have length 1.")
  }
  
  Y <- data[data>t]-t
  n <- length(Y)
  
  if (abs(gamma) < .Machine$double.eps) {
    R <- Y/sigma
  } else {
    R <- 1/gamma*log(1+gamma/sigma*Y)
  }
  
  # calculate theoretical and empirical quantiles)
  res.the <- -log( 1 - (1:n) / (n+1))
  res.emp<- sort(R)
  
  
  # plots if TRUE
  .plotfun(res.the, res.emp, type="p", xlab="Quantiles of Standard Exponential", ylab="R", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles gqq.the
  # and empirical quantiles gqq.emp
  .output(list(res.the=res.the, res.emp=res.emp), plot=plot, add=FALSE)
  
  
}
