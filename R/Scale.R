
#Scale estimator
Scale <- function(data, gamma = NULL, plot = FALSE, add = FALSE, 
                    main = "Estimates of scale parameter", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  if(is.null(gamma)) {
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
  .plotfun(K, A[K], type="l", xlab="k", ylab="scale", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K,A=A[K],C=C[K]), plot=plot, add=add)
  
}


#Bias-reduced scale estimator
Scale.2o <- function(data, gamma, b, beta, plot = FALSE, add = FALSE,
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
  .plotfun(K, A[K], type="l", xlab="k", ylab="scale", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K,A=A[K],C=C[K]), plot=plot, add=add)
}