
# This file contains the implementation of the generalised Hill estimator.

###########################################################################

# Computes the Hill estimates of gamma (Section 5.2.3 in Beirlant et al. (2004)) 
# for a numeric vector of observations (data) and as a
# function of k
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot


genHill <- function(data, gamma, logk = FALSE, plot = FALSE, add = FALSE,  
                    main = "Generalised Hill estimates of the EVI", ...) {
  
  # Check input arguments
  .checkInput(data, gamma)
  
  n <- length(data)
  UH.scores <- numeric(n)
  Hill <- numeric(n)
  ind <- 1:(n-1)
  X <- as.numeric(sort(data))
  
  
  ### Hill estimates
  
  ######################
  # Vectorised version
  UH.scores[ind] <- X[n-ind] * gamma[ind]
  
  # Stop at n-2 since UH.scores[n]=0, so log(UH.scores[n])=Inf
  K <- 1:max((n-2), 1) 
  Hill[K] <- cumsum(log(UH.scores[K])) / K - log(UH.scores[K+1])
  
  ######################
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), Hill[K], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }
 
  # output list with values of k and corresponding Hill estimates
  .output(list(k=K, gamma=Hill[K]), plot=plot, add=add)
  
}


#Estimate extreme quantiles using the generalised Hill estimator for gamma
#p. 156-157 in Beirlant et al. (2004)
QuantGH <- function(data, gamma, p, plot = FALSE, add = FALSE, 
                    main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data, gamma, gammapos=FALSE)
  .checkProb(p)
  
  
  X <- sort(data)
  n <- length(X)
  a <- numeric(n)
  quant <- numeric(n)
  K <- 1:length(gamma)
  
  # Estimator for extreme quantiles
  
  H <- Hill(X)$gamma
  a <- X[n-K]*H[K]*(1-pmin(gamma, 0))
  
  quant[K] <- X[n-K] + a/gamma[K]*( ((K+1)/((n+1)*p))^gamma[K] - 1 )
  
  # plots if TRUE
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p), plot=plot, add=add)
  
}


ProbGH <- function(data, gamma, q, plot = FALSE, add = FALSE, 
                   main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data, gamma, gammapos=FALSE)
  
  if (length(q) > 1) {
    stop("q should be a numeric of length 1.")
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  prob <- numeric(n)
  K <- 1:length(gamma)
  
  # Estimator for tail probabilities
  
  H <- Hill(X)$gamma
  a <- X[n-K]*H[K]*(1-pmin(gamma, 0))
  
  prob[K] <- ((K+1)/(n+1)) * (1 + gamma[K]/a[K]*(q-X[n-K]))^(-1/gamma[K])
  prob[prob < 0 | prob > 1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding probability estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}


ReturnGH <- function(data, gamma, q, plot = FALSE, add = FALSE, 
                     main = "Estimates of large return period", ...) {
  
  # Check input arguments
  .checkInput(data, gamma, gammapos=FALSE)
  
  if (length(q) > 1) {
    stop("q should be a numeric of length 1.")
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  r <- numeric(n)
  K <- 1:length(gamma)
  
  # Estimator for tail probabilities
  
  H <- Hill(X)$gamma
  a <- X[n-K]*H[K]*(1-pmin(gamma, 0))
  
  r[K] <- (n+1)/(K+1) * (1 + gamma[K]/a[K]*(q-X[n-K]))^(1/gamma[K])
  r[r < 0] <- NA
  
  # plots if TRUE
  .plotfun(K, r[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=r[K], q=q), plot=plot, add=add)
  
}
