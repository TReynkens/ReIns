
# Kaplan-Meier estimator for the CDF evaluated in x
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary sample CDF
#
# Pure R implementation
.KaplanMeier_R <- function(x, data, censored) {
  
  
  # Check input arguments
  if (!is.numeric(data)) {
    stop("data should be a numeric vector.")
  }
  if (!is.numeric(x)) {
    stop("x should be a numeric.")
  }
  censored <- .checkCensored(censored, length(data))
  
  
  # Sort the data with index return
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  
  est <- numeric(length(x))
  
  for (i in 1:length(x)) {
    est[i] <- prod(1-delta[X<=x[i]]/(n-(1:n)[X<=x[i]]+1))
  }
  
  # Kaplan Meier estimator for the CDF
  return( 1 - est )
}


# Kaplan-Meier estimator for the CDF evaluated in x
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary sample CDF
#
.KaplanMeier_surv <- function(x, data, censored, conf.type = "plain") {
  
  
  # Check input arguments
  if (!is.numeric(data)) {
    stop("data should be a numeric vector.")
  }
  if (!is.numeric(x)) {
    stop("x should be a numeric.")
  }
  censored <- .checkCensored(censored, length(data))
  
  
  # Sort the data with index return
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  
  fit  <- survfit(Surv(X, delta) ~1, conf.type=conf.type)

  f = stepfun(fit$time,c(1,fit$surv))
  est <- f(x)
  
  
  # Kaplan Meier estimator for the CDF
  L <- list(km = 1 - est)

  # Lower CI
  if(!is.null(fit$lower)) {
    fl <- stepfun(fit$time,c(1,fit$lower))
    clower <- fl(x)
    L$clower <- 1-clower
  }


  # Upper CI
  if(!is.null(fit$upper)) {
    fu <- stepfun(fit$time,c(1,fit$upper))
    cupper <- fu(x)
    L$cupper <- 1-cupper
  } 
  
 return(L)
}




# Use version based on survival package by default when x is a large vector
KaplanMeier <- function(x, data, censored) {
  
 if (length(x)<=1000) {
   return(.KaplanMeier_R(x, data, censored))
 } else {
   return(.KaplanMeier_surv(x, data, censored)$km)
 }
}


# Turnbull estimator for the CDF evaluated in x
# L and R are the lower and upper values for interval censoring
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary sample CDF
#
Turnbull <- function(x, L, R, censored, trunclower = 0, truncupper = Inf, conf.type = "plain", conf.int = 0.95) {
  
  
  # Check input arguments
  if (!is.numeric(L)) {
    stop("L should be a numeric vector.")
  }
  
  # Check input arguments
  if (!is.numeric(R)) {
    stop("R should be a numeric vector.")
  }
  
  if (!is.numeric(trunclower)) {
    stop("trunclower should be numeric.")
  }
  
  if (!is.numeric(truncupper)) {
    stop("truncupper should be numeric.")
  }
  
  # Check lengths
  if (length(L)!=length(R)) {
    stop("L and R should have the same length.")
  }
  
  if (length(trunclower)!=1) {
    stop("trunclower should have length 1.")
  }
  
  if (length(truncupper)!=1) {
    stop("truncupper should have length 1.")
  }
  
  if (any(L<trunclower)) {
    stop("All elements of L should be larger than trunclower.")
  }
  
  if (any(R>truncupper)) {
    stop("All elements of R should be smaller than truncupper.")
  }
  
  if (any(L>R)) {
    stop("Each element of L should be smaller or equal than each corresponding element of R.")
  }
  
  if (!is.numeric(x)) {
    stop("x should be numeric.")
  }
  censored <- .checkCensored(censored, length(L))
  
  
  ft <- .Turnbull_internal(L=L, R=R, censored=censored, conf.type = conf.type, conf.int = conf.int,
                           trunclower=trunclower, truncupper=truncupper)
  est <- ft$f(x)
  
  # Turnbull estimator for the CDF
  return( list(cdf = 1 - est, fit=ft$fit) )
}


# Function to return Turnbull stepfunction for survival function
.Turnbull_internal <- function(L, R, censored, trunclower = 0, truncupper = Inf, 
                               conf.type = "plain", conf.int = 0.95) {
  
  # Uncensored observations
  R[!censored] <- L[!censored]
  # Right censored observations, set upper bound to NA
  R[censored & R==truncupper] <- NA
  # Left censored observations, set lower bound to NA
  L[censored & L==trunclower] <- NA


  fit  <- survfit(Surv(time=L, time2=R, type="interval2") ~1,
                  conf.type = conf.type, conf.int = conf.int)

  
  f = stepfun(fit$time,c(1,fit$surv))
              
  return(list(f=f, fit=fit))
}