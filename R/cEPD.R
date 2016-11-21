
# This file contains the implementation of the EPD estimator adapted for right censored data.

###########################################################################


# Proportion of non-censored observations
.p_hat_fun <- function(delta.n) {
  
  n <- length(delta.n)
  K <- 1:(n-1)
  
  p.hat <- cumsum(delta.n[n-K+1])/K
#   for (k in 1:(n-1)) {
#     p.hat[k] <- mean(delta.n[n:(n-k+1)])
#   }
  return(p.hat)
}



# E_{Z,k}(s)
.ES <- function(s, Z.n) {
  
  n <- length(Z.n)
  output <- numeric(n-1)
  
  for (k in 1:(n-1)) {
    v <- (Z.n[n-(1:k)+1]/Z.n[n-k])^(s[k])
    output[k] <- mean(v)
  }
  return(output)
}



# E^(c)_{Z,k}(s)
.cEs <- function(s, Z.n, delta.n) {
  
  n <- length(Z.n)
  output <- numeric(n-1)
  
  for (k in 1:(n-1)) {
    v <- delta.n[n-(1:k)+1] * (Z.n[n-(1:k)+1]/Z.n[n-k])^(s[k])
    output[k] <- mean(v)
  }
  return(output)
}



# Censored EPD estimator for gamma_1 and kappa_1
#
# rho is a parameter for the rho_estimator of Fraga Alves et al. (2003)
# when strictly positive or (a) choice(s) for rho if negative

cEPD <- function(data, censored, rho = -1, beta = NULL, logk = FALSE, plot = FALSE, add = FALSE, 
                 main = "EPD estimates of EVI", ...) {
  
  # Check input
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  n <- length(data)
  k <- n-1 
  
  s <- sort(data, index.return = TRUE)
  Z.n <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta.n <- !(censored[s$ix])
  
  HillZ <- Hill(Z.n)$gamma
  phat <- .p_hat_fun(delta.n)
  
  nrho <- length(rho)

  if (is.null(beta)) {
    # determine beta using rho
    
    beta <- matrix(0, n-1, nrho)
    
    if (all(rho > 0) & nrho == 1) {
      rho <- .rhoEst(data, alpha=1, tau=rho)$rho
     
      # Estimates for rho of Fraga Alves et al. (2003) used 
      # and hence a different value of beta for each k
      beta[,1] <- -rho/HillZ
      
      for(j in 1:length(rho)) {
        beta[,j] <- -rho[j]/HillZ
      }
      
    } else if (all(rho < 0)) {
      
      # rho is provided => beta is constant over k
      # (but differs with rho)
      for(j in 1:nrho) {
        beta[,j] <- -rho[j]/HillZ
      }
      
      
    } else {
      stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }

  } else {
    nrho <- length(beta)
    
    # use provided value(s) for beta
    if (length(beta) == 1) {
      beta <- matrix(beta, n-1, nrho)
#     } else if (length(beta) == n-1) {
#       beta <- as.matrix(rep(beta,nrho),ncol=length(rho),byrow=FALSE)
#     } else {
#       stop(paste0("beta should have length 1 or n-1 = ",n-1,"."))
#     }
    } else {
      beta <- matrix(rep(beta, n-1), ncol=length(beta), byrow=TRUE)
    }
  }
  
  gamma1 <- matrix(0, n-1, nrho)
  kappa1 <- matrix(0, n-1, nrho)
  Delta <- matrix(0, n-1, nrho)
  
  for(j in 1:nrho) {

    K <- 1:k
    D <- - (beta[K,j]^4 * HillZ[K]^3) / ( (1+HillZ[K]*beta[K,j])^2 * (1+2*HillZ[K]*beta[K,j]) )
    
    Es <- .ES(-beta[,j], Z.n)
    cEs <- .cEs(-beta[,j], Z.n, delta.n)

    # Estimates for kappa1
    kappa1[,j] <- (1 - Es[K] - beta[K,j] * (HillZ[K] / phat[K]) * cEs[K]) / D[K]
    kappa1[,j] <- pmax(kappa1[,j], pmax(-1, -1/beta[,j])+0.001)

    # Estimates for Delta
    Delta[,j] <- (kappa1[,j]*(1-Es))/phat 
    
    # Has to be strictly negative otherwise the bias is increased,
    # hence we put it to 0 such that gamma1 = cHill then
    Delta[,j] <- pmin(0, Delta[,j])
    
    # gamma1 = cHill + Delta
    gamma1[,j] <- (HillZ/phat) + Delta[,j]
    gamma1[gamma1[,j] <= 0, j] <- 0.001

  }
  # Only positive values are allowed
#   gamma1[gamma1 <= 0] <- NA
#   
# 
#   kappa1[kappa1 <= pmax(-1,-1/beta)] <- NA
  
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), gamma1[,1], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, gamma1[,1], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }
  
  # Transform to vectors if rho is a single value
  if (nrho == 1) {
    gamma1 <- as.vector(gamma1)
    kappa1 <- as.vector(kappa1)
    beta <- as.vector(beta)
    Delta <- as.vector(Delta)
  } else if (plot | add) {
  # Add lines
    for(j in 2:nrho) {
      lines(K, gamma1[,j], lty=j)
    }
  }
  
  .output(list(k=K, gamma1=gamma1, kappa1=kappa1, beta=beta, Delta=Delta), plot=plot, add=add)
}



# Estimator for small exceedance probabilities for (right) censored data
# using EPD estimates for (right) censored data
cProbEPD <- function(data, censored, gamma1, kappa1, beta, q, plot = FALSE, add=FALSE,
                    main = "Estimates of small exceedance probability", ...) {

  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q) > 1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- as.numeric(s$x)
  sortix <- s$ix
  n <- length(X)
  prob <- rep(NA, n)
  K <- 1:(n-1)
  
  K2 <- K[!is.na(gamma1[K])]
  
  # Kaplan-Meier estimator for CDF in X[n-K]
  km <- KaplanMeier(X[n-K2], data=X, censored = censored[sortix])$surv

  prob[K2] <- km * (1-pepd(q/X[n-K2], gamma=gamma1[K2], kappa=kappa1[K2], tau=-beta[K2]))
  prob[prob < 0 | prob > 1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}


# Estimator for return periods for (right) censored data
# using EPD estimates for (right) censored data
cReturnEPD <- function(data, censored, gamma1, kappa1, beta, q, plot = FALSE, add = FALSE,
                     main = "Estimates of large return period", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q) > 1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- as.numeric(s$x)
  sortix <- s$ix
  n <- length(X)
  R <- rep(NA, n)
  K <- 1:(n-1)
  
  K2 <- K[!is.na(gamma1[K])]
  
  # Kaplan-Meier estimator for CDF in X[n-K]
  km <- KaplanMeier(X[n-K2], data=X, censored = censored[sortix])$surv
  
  R[K2] <- 1 / (km * (1-pepd(q/X[n-K2], gamma=gamma1[K2], kappa=kappa1[K2], tau=-beta[K2])))
  R[R < 1] <- NA
  
  # plots if TRUE
  .plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=R[K], q=q), plot=plot, add=add)
  
}
