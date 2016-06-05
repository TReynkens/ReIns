

# Regression on scale parameter with constant gamma (estimator for A(s))
#
# s is the point in which the estimated scale function is evaluated
#
# Z is the response variable data
#
# kernel is the kernel function used and h is the bandwidth
# If h is equal to NULL, a bandwidth will be selected (similar to density).

ScaleReg <- function(s, Z, kernel =  c("gaussian", "epanechnikov", "rectangular",
                                       "triangular", "biweight", "cosine", "optcosine"), 
                     h = NULL, plot = TRUE, add = FALSE, 
                     main = "Estimates of scale parameter", ...) {
  
  # Check input arguments
  .checkInput(Z)
  
  # Sort data
  Zsort <- as.numeric(sort(Z))
  n <- length(Z)
  K <- 1:(n-1)
  
  # Select right kernel function with bandwidth h (K_h(x) = K(x/h)/h)
  kernel <- match.arg(kernel)
  kernelh <- .kernel_aux(kernel=kernel, h=h)
  
  # Scale estimates in the point s for several values of k
  A <- numeric(n)
  
  for(k in K) {
    A[k] <- 1/(k+1) * sum( (Z>Zsort[n-k]) * kernelh( s-(1:n)/n ) )
  }
  
  # plots if TRUE
  .plotfun(K, A[K], type="l", xlab="k", ylab="Scale", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K,A=A[K]), plot=plot, add=add)
}

# Estimator of small tail probability 1-F_i(x) in regression case where gamma is constant
# A are the estimates for A(i/n) for k = 1, ...,n-1
ProbReg <- function(Z, A, q, plot = FALSE, add = FALSE, 
                    main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(Z, scale=A)
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  
  Zsort <- as.numeric(sort(Z))
  n <- length(Z)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Hill estimator for gamma
  hill <- Hill(Z, plot=FALSE)$gamma
  
  # Estimator for small exceedance probability
  prob[K] <- A[K] * (K+1)/(n+1) * (q/Zsort[n-K])^(-1/hill[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
}


# Estimator of large quantile Q_i(1-p) in regression case where gamma is constant
# A are the estimates for A(i/n) for k = 1, ...,n-1
QuantReg <- function(Z, A, p, plot = FALSE, add = FALSE, 
                     main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(Z, scale=A)
  
  .checkProb(p)
  
  Zsort <- as.numeric(sort(Z))
  n <- length(Z)
  Q <- numeric(n)
  K <- 1:(n-1)
  
  # Hill estimator for gamma
  hill <- Hill(Z, plot=FALSE)$gamma
  
  # Estimator for quantiles
  Q[K] <- Zsort[n-K] * ((K+1)/((n+1)*p) * A[K])^(hill[K])
  
  # plots if TRUE
  .plotfun(K, Q[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=Q[K], p=p),plot=plot,add=add)
}


# Auxiliary function to select kernel and returns K_h(x) = kernel(x/h)/h
# h is the chosen bandwidth and if equal to NULL, a bandwidth will be selected
# (similar to density)
.kernel_aux <- function(kernel =  c("gaussian", "epanechnikov", "rectangular",
                                   "triangular", "biweight", "cosine", "optcosine"), h = NULL) {
  
  kernel <- match.arg(kernel)
  
  # Set default bandwidth
  if(is.null(h) | !is.numeric(h)) {
    bw <- switch(kernel,
                 gaussian = 1/(2*sqrt(pi)),
                 rectangular = sqrt(3)/6,
                 triangular  = sqrt(6)/9,
                 epanechnikov= 3/(5*sqrt(5)),
                 biweight    = 5*sqrt(7)/49,
                 cosine      = 3/4*sqrt(1/3 - 2/pi^2),
                 optcosine   = sqrt(1-8/pi^2)*pi^2/16
    )
  } else {
    bw <- h
  }
  
  # Select kernel function
  kernelfun <- switch(kernel,
                  gaussian = function(x) { dnorm(x, sd = bw) },
                  ## In the following, a := bw / sigma(K0), where
                  ##	K0() is the unscaled kernel below
                  rectangular =  function(x) {
                    a <- bw*sqrt(3)
                    return( ifelse(abs(x) < a, .5/a, 0)  ) },
                  triangular =  function(x) {
                    a <- bw*sqrt(6) ; ax <- abs(x)
                    return( ifelse(ax < a, (1 - ax/a)/a, 0)  ) },
                  epanechnikov =  function(x) {
                    a <- bw*sqrt(5) ; ax <- abs(x)
                    return( ifelse(ax < a, 3/4*(1 - (ax/a)^2)/a, 0)  ) },
                  biweight =  function(x) { ## aka quartic
                    a <- bw*sqrt(7) ; ax <- abs(x)
                    return( ifelse(ax < a, 15/16*(1 - (ax/a)^2)^2/a, 0)  ) },
                  cosine =  function(x) {
                    a <- bw/sqrt(1/3 - 2/pi^2)
                    return( ifelse(abs(x) < a, (1+cos(pi*x/a))/(2*a),0) ) },
                  optcosine =  function(x) {
                    a <- bw/sqrt(1-8/pi^2)
                    return( ifelse(abs(x) < a, pi/4*cos(pi*x/(2*a))/a, 0) ) }
  )
  
  # K_h(x) = K(x/h)/h
  kh <- function(x) { kernelfun(x/bw) / bw}
  
  return(kh)
}