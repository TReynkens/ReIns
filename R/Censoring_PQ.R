
# This file contains the estimators for quantiles, small probabilities and return periods
# suitable for right censored data.

###########################################################################


# Computes estimates of small exceedance probability 1-F(q) 
# for a numeric vector of observations 
# (data) and as a function of k
#
# Estimates are based on prior estimates gamma1 for the
# extreme value index (cHill)
#
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results indicates that non of the data points is censored
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot

cProb <- function(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
                        main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data,gamma1)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  #delta <- !(censored[sortix]*1)
  n <- length(X)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Kaplan-Meier estimator for survival function in X[n-K]
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  # Weissman estimator for probabilities
  prob[K] <- km * (q/X[n-K])^(-1/gamma1[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
}


cReturn <- function(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
                  main = "Estimates of large return period", ...) {
  
  # Check input arguments
  .checkInput(data,gamma1)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  #delta <- !(censored[sortix]*1)
  n <- length(X)
  R <- numeric(n)
  K <- 1:(n-1)
  
  # Kaplan-Meier estimator for survival function in X[n-K]
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  # Weissman estimator for probabilities
  R[K] <- 1 / ( km * (q/X[n-K])^(-1/gamma1[K]) )
  R[R<1] <- NA
  
  # plots if TRUE
  .plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=R[K], q=q),plot=plot,add=add)
}



##########################################################################

# Computes estimates of extreme quantile Q(1-p)  for a numeric vector of observations 
# (data) and as a function of k
#
# Estimates are based on prior estimates gamma1 for the
# extreme value index (cHill)
#
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results indicates that non of the data points is censored
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot

cQuant <- function(data, censored, gamma1, p, plot = FALSE, add = FALSE, 
                        main = "Estimates of extreme quantile", ...) {
  
  
  # Check input arguments
  .checkInput(data,gamma1)
  censored <- .checkCensored(censored, length(data))
  .checkProb(p)
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  #delta <- !(censored[sortix]*1)
  n <- length(X)
  quant <- numeric(n)
  K <- 1:(n-1)
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  # Weisman estimator for quantiles
  quant[K] <- X[n-K] * (km/p)^(gamma1[K])
  
  # plots if TRUE
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p),plot=plot,add=add)
}


#################################################################################

#Estimate extreme quantiles using the censored generalised Hill estimator for gamma1
cQuantGH <- function(data, censored, gamma1, p, plot = FALSE, add = FALSE, 
                     main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data,gamma1,gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  .checkProb(p)
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  quant <- numeric(n)
  K <- 1:length(gamma1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Hill estimator (non-censored)
  H <- Hill(X)$gamma
  
  # Auxiliary value
  a <- X[n-K] * H[K] * (1-pmin(gamma1[K],0)) / pk[K]
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  
  # Estimator for extreme quantiles
  quant[K] <- X[n-K] + a/gamma1[K] * ( (km/p)^gamma1[K] - 1 )
  
  # plots if TRUE
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p), plot=plot, add=add)
  
}



#Estimate small exceedance probabilities using the censored generalised Hill estimator for gamma1
cProbGH <- function(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
                    main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data,gamma1,gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  prob <- numeric(n)
  K <- 1:length(gamma1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Hill estimator (non-censored)
  H <- Hill(X)$gamma
  
  # Auxiliary value
  a <- X[n-K] * H[K] * (1-pmin(gamma1[K],0)) / pk[K]
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  prob[K] <- km * (1 + gamma1[K]/a[K]*(q-X[n-K]))^(-1/gamma1[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding probability estimates
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}


#Estimate return periods using the censored generalised Hill estimator for gamma1
cReturnGH <- function(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
                    main = "Estimates of large return period", ...) {
  
  # Check input arguments
  .checkInput(data,gamma1,gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  R <- numeric(n)
  K <- 1:length(gamma1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Hill estimator (non-censored)
  H <- Hill(X)$gamma
  
  # Auxiliary value
  a <- X[n-K] * H[K] * (1-pmin(gamma1[K],0)) / pk[K]
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  R[K] <- 1 / ( km * (1 + gamma1[K]/a[K]*(q-X[n-K]))^(-1/gamma1[K]) )
  R[R<1] <- NA
  
  # plots if TRUE
  .plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=R[K], q=q), plot=plot, add=add)
  
}


#Same expressions using the moment estimator
cQuantMOM <- cQuantGH
cProbMOM <- cProbGH
cReturnMOM <- cReturnGH


#Estimate extreme quantiles using the censored POT estimator for gamma1
cQuantGPD <- function(data, censored, gamma1, sigma1, p, plot = FALSE, add = FALSE, 
                      main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data, gamma1, scale=sigma1, gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  .checkProb(p)
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  quant <- numeric(n)
  K <- 1:(n-1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Auxiliary value
  a <- sigma1[K]/pk
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  
  # Estimator for extreme quantiles
  quant[K] <- X[n-K] + a/gamma1[K] * ( (km/p)^gamma1[K] - 1 )
  
  # plots if TRUE
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p), plot=plot, add=add)
  
}

#Estimate small exceedance probabilities using the censored POT estimator for gamma1
cProbGPD <- function(data, censored, gamma1, sigma1, q, plot = FALSE, add = FALSE, 
                     main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data, gamma1, scale=sigma1, gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  prob <- numeric(n)
  K <- 1:(n-1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Auxiliary value
  a <- sigma1[K]/pk
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  

  prob[K] <- km * (1 + gamma1[K]/a[K]*(q-X[n-K]))^(-1/gamma1[K])
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding probabilities
  # and the considered large quantile q
  
  .output(list(k=K, P=prob[K], q=q), plot=plot, add=add)
  
}


#Estimate return period using the censored POT estimator for gamma1
cReturnGPD <- function(data, censored, gamma1, sigma1, q, plot = FALSE, add = FALSE, 
                     main = "Estimates of large return period", ...) {
  
  # Check input arguments
  .checkInput(data, gamma1, scale=sigma1, gammapos=FALSE)
  censored <- .checkCensored(censored, length(data))
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  R <- numeric(n)
  K <- 1:(n-1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Auxiliary value
  a <- sigma1[K]/pk
  
  # Kaplan-Meier estimator for CDF
  km <- KaplanMeier(X[n-K], data = X, censored = censored[sortix])$surv
  
  
  R[K] <- 1 / ( km * (1 + gamma1[K]/a[K]*(q-X[n-K]))^(-1/gamma1[K]) )
  R[R<1] <- NA
  
  # plots if TRUE
  .plotfun(K, R[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=R[K], q=q), plot=plot, add=add)
  
}