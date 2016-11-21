
# This file contains implementations of scale estimators.

###########################################################################

# Scale estimator
Scale <- function(data, gamma = NULL, logk = FALSE, plot = FALSE, add = FALSE, 
                    main = "Estimates of scale parameter", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  if (is.null(gamma)) {
    gamma <- Hill(data)$gamma
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  
  A <- numeric(n)
  C <- numeric(n)
  
  A[K] <- (K+1)/(n+1) * X[n-K]^(1/gamma[K])
  C[K] <- A[K]^gamma[K]
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), A[K], type="l", xlab="log(k)", ylab="scale", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, A[K], type="l", xlab="k", ylab="scale", main=main, plot=plot, add=add, ...)
  }
  
  .output(list(k=K, A=A[K], C=C[K]), plot=plot, add=add)
  
}


# Bias-reduced scale estimator using Hill.2oQV
Scale.2o <- function(data, gamma, b, beta, logk = FALSE, plot = FALSE, add = FALSE,
                       main = "Estimates of scale parameter", ...) {
  
  # Check input arguments
  .checkInput(data)

  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  A <- numeric(n)
  C <- numeric(n)
  
  #beta=-rho here
  A[K] <- (K+1)/(n+1) * X[n-K]^(1/gamma[K]) * exp(b[K]/beta[K])
  C[K] <- A[K]^gamma[K]
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), A[K], type="l", xlab="log(k)", ylab="scale", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, A[K], type="l", xlab="k", ylab="scale", main=main, plot=plot, add=add, ...)
  }
  
  .output(list(k=K, A=A[K], C=C[K]), plot=plot, add=add)
}


# Bias-reduced scale estimator using EPD
ScaleEPD <- function (data, gamma, kappa, logk = FALSE, plot = FALSE, add = FALSE, 
                      main = "Estimates of scale parameter", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  K <- 1:(n-1)
  A <- numeric(n)
  C <- numeric(n)
  
  #   kappa[kappa<pmax(-1,-1/gamma[K])] <- NA
  gamma[gamma <= 0] <- NA
  
  A[K] <- (K + 1)/(n + 1) * X[n - K]^(1/gamma[K]) * (1-kappa[K]/gamma[K])
  C[K] <- A[K]^gamma[K]
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), A[K], type="l", xlab="log(k)", ylab="scale", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, A[K], type="l", xlab="k", ylab="scale", main=main, plot=plot, add=add, ...)
  }
  
  .output(list(k=K, A=A[K], C=C[K]), plot=plot, add=add)
}

