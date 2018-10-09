
# This file contains EVT methods suitable for right censored data.
# Estimators for quantiles, small probabilities and return periods are in Censoring_PQ.R.
# Estimators suitable for interval censored data are in IntervalCensoring.R

#####################################################################

# Exponential QQplot for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 indicates that non of the data points is censored
cExpQQ <- function(data, censored, plot = TRUE, main = "Exponential QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  n <- length(X)
  
  
  K <- 1:(n-1)
  
  # calculate theoretical and empirical quantiles
  
  
  # -log of Kaplan-Meier estimator for the survival function in Z_[n-K+1]
  eqq.the <- -log(KaplanMeier(X[n-K+1], data = X, censored = censored[sortix])$surv )
  eqq.emp <- X[n-K+1]
  
  
  # plots if TRUE
  .plotfun(eqq.the, eqq.emp, type="p", xlab="-log of Kaplan-Meier estimator", ylab="Z", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles eqq.the and empirical quantiles eqq.emp
  .output(list(eqq.the=eqq.the, eqq.emp=eqq.emp), plot=plot, add=FALSE)
  
}

#Log-normal QQplot
cLognormalQQ <- function(data, censored, plot = TRUE, main = "Log-normal QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  
  
  # calculate theoretical and empirical quantiles
  
  #  qnorm of Kaplan-Meier estimator for the survival function in Z_[i]
  lnqq.the <- qnorm(1-KaplanMeier(X, data = X, censored = censored[sortix])$surv )
  lnqq.emp <- log(X)
  
  # plots if TRUE
  .plotfun(lnqq.the, lnqq.emp, type="p", xlab="Normal quantile of Kaplan-Meier estimator", ylab="log(Z)", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles lnqq.the and empirical quantiles lnqq.emp
  .output(list(lnqq.the=lnqq.the, lnqq.emp=lnqq.emp), plot=plot, add=FALSE)
  
}

# Pareto QQplot for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 indicates that non of the data points is censored
cParetoQQ <- function(data, censored, plot = TRUE, main = "Pareto QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  n <- length(X)
  
  
  K <- 1:(n-1)
  
  # calculate theoretical and empirical quantiles
  
  
  # -log of Kaplan-Meier estimator for the survival function in Z_[n-K+1]
  pqq.the <- -log( KaplanMeier(X[n-K+1], data = X, censored = censored[sortix])$surv )
  pqq.emp <- log(X[n-K+1])
  
  
  # plots if TRUE
  .plotfun(pqq.the, pqq.emp, type="p", xlab="-log of Kaplan-Meier estimator", ylab="log(Z)", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(pqq.the=pqq.the, pqq.emp=pqq.emp), plot=plot, add=FALSE)
  
}

# Weibull QQ-plot
cWeibullQQ <- function(data, censored, plot = TRUE, main = "Weibull QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  n <- length(X)
  
  
  K <- 1:(n-1)
  
  # calculate theoretical and empirical quantiles
  
  # log of -log of Kaplan-Meier estimator for the survival function in Z_[n-K+1]
  wqq.the <- log(-log(KaplanMeier(X[n-K+1], data = X, censored = censored[sortix])$surv))
  wqq.emp <- log(X[n-K+1])
  
  # plots if TRUE
  .plotfun(wqq.the, wqq.emp, type="p", xlab="log of -log of Kaplan-Meier estimator", ylab="log(Z)", 
           main=main, plot=plot, add=FALSE, ...)
  
  # output list with theoretical quantiles wqq.the and empirical quantiles wqq.emp
  .output(list(wqq.the=wqq.the, wqq.emp=wqq.emp), plot=plot, add=FALSE)
  
}


##################################################################################


# Hill estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary Hill estimator
cHill <- function(data, censored, logk = FALSE, plot = FALSE, add = FALSE, main = "Hill estimates of the EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  if (length(censored) != 1) {
    if (length(data) != length(censored)) {
      stop("data and censored should have the same length.")
    }
  } else {
    censored <- rep(0, length(data))
  }
  
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[s$ix])
  n <- length(X)
  Hill <- numeric(n)
  
  if (n == 1) {
    stop("We need at least two data points.")
  }
  
  K <- 1:(n-1)
  
  # Hill estimates
  # Fast vectorised version
  Hill[K] <- (cumsum(log(X[n-K+1]))- K*log(X[n-K])) / cumsum(delta[n-K+1]) 
  # Problems with division by 0
  Hill[K[cumsum(delta[n-K+1]) == 0]] <- NA
  

  # plots if TRUE
  if (logk) {
    .plotfun(log(K), Hill[K], type="l", xlab="log(k)", ylab="gamma1", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  }
  
  # output list with values of k and corresponding Hill estimates
  
  .output(list(k=K, gamma1=Hill[K]), plot=plot, add=add)
  
}


# generalised Hill estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary generalised Hill estimator.
cgenHill <- function(data, censored, logk = FALSE, plot = FALSE, add = FALSE,  
         main = "Generalised Hill estimates of the EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  UH.scores <- numeric(n)
  Hill <- numeric(n)
  ind <- 1:(n-1)
  
  if (n == 1) {
    stop("We need at least two data points.")
  }
  
  ### Generalised Hill estimates
  
  # Fast vectorised version
  UH.scores[ind] <- X[n-ind] * Hill(X)$gamma[ind]

  # Stop at n-2 since UH.scores[n]=0, so log(UH.scores[n])=Inf
  K <- 1:max((n-2), 1) 
  Hill[K] <- (cumsum(log(UH.scores[K]))- K*log(UH.scores[K+1])) / cumsum(delta[n-K+1]) 
  # Problems with division by 0
  Hill[K[cumsum(delta[n-K+1]) == 0]] <- NA
  

  # plots if TRUE
  if (logk) {
    .plotfun(log(K), Hill[K], type="l", xlab="log(k)", ylab="gamma1", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  }
  
  # output list with values of k and
  # corresponding generalised Hill estimates
  
  .output(list(k=K, gamma1=Hill[K]), plot=plot, add=add)
  
}

# Moment estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary Moment estimator
cMoment <- function(data, censored, logk = FALSE, plot = FALSE, add = FALSE, main = "Moment estimates of the EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[s$ix])
  n <- length(X)
  M1 <- numeric(n)
  M2 <- numeric(n)
  Mom <- numeric(n)
  K <- 1:(n-1)
  
  
  
  ### Moment estimates
  
  ######################
  #Fast vectorised version
  
  M1[K] <- cumsum(log(X[n-K+1]))/K - log(X[n-K])
  # for k=1:(n-1)
  # 1/k*sum_{j=1}^k [log(X[n-j+1])-log(X[n-k])]^2
  # = 1/k*sum_{j=1}^k log(X[n-j+1])^2 - 2*log(X[n-k])*1/k*sum_{j=1}^k log(X[n-j+1]) + log(X[n-k])^2
  M2[K] <- cumsum(log(X[n-K+1])^2)/K-2*cumsum(log(X[n-K+1]))/K*log(X[n-K]) + log(X[n-K])^2
  
  Mom <- M1 + 1 - (1-M1^2/M2)^(-1)/2
  #Add censoring
  Mom[K] <- Mom[K] / (cumsum(delta[n-K+1])/K )

  # Problems with division by 0
  Mom[K[cumsum(delta[n-K+1]) == 0]] <- NA
  
  ######################
  
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), Mom[K], type="l", xlab="log(k)", ylab="gamma1", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, Mom[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  }
 
  # output list with values of k and
  # corresponding estimates for gamma1, b and beta
  .output(list(k=K, gamma1=Mom[K]), plot=plot, add=add)
  
}


cGPDmle <- function(data, censored, start = c(0.1, 1), warnings = FALSE, logk = FALSE, 
                 plot = FALSE, add = FALSE, main = "POT estimates of the EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  censored <- .checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  K <- 1:(n-1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Apply ordinary POT estimator to the data
  L <- POT(data=X, start=start, warnings=warnings, plot=FALSE, add=FALSE)
  
  
  # Estimates
  gamma1 <- L$gamma[K]
  sigma1 <- L$sigma[K]
  
  # Censored version
  gamma1 <- gamma1 / pk
    
  # plots if TRUE
  if (logk) {
    .plotfun(log(K), gamma1[K], type="l", xlab="log(k)", ylab="gamma1", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, gamma1[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  }
 
  
  .output(list(k=K, gamma1=gamma1[K], sigma1=sigma1[K]), plot=plot, add=add)
}

cPOT <- cGPDmle

