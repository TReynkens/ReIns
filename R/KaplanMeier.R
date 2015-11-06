
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
Turnbull <- function(x, L, R, censored, conf.type = "plain", conf.int = 0.95) {
  
  
  # Check input arguments
  if (!is.numeric(L)) {
    stop("L should be a numeric vector.")
  }
  
  # Check input arguments
  if (!is.numeric(R)) {
    stop("R should be a numeric vector.")
  }
  
  
  # Check input arguments
  if (length(L)!=length(R)) {
    stop("L and R should have the same length.")
  }
  
  
  if (!is.numeric(x)) {
    stop("x should be numeric.")
  }
  censored <- .checkCensored(censored, length(L))
  
  
  ft <- .Turnbull_internal(L=L, R=R, censored=censored, conf.type = conf.type, conf.int = conf.int)
  est <- ft$f(x)
  
  # Turnbull estimator for the CDF
  return( list(cdf = 1 - est, fit=ft$fit) )
}


# Function to return Turnbull stepfunction for survival function
# No left censoring included!
.Turnbull_internal <- function(L, R, censored, conf.type = "plain", conf.int = 0.95) {
  
  event <- !censored
  # # 3 for interval censoring and 0 for right censoring
  event[event==0 & is.finite(R)] <- 3
  
  type <- "interval"

  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, 
                  conf.type = conf.type, conf.int = conf.int)

  
  f = stepfun(fit$time,c(1,fit$surv))
              
  return(list(f=f, fit=fit))
}


# # Sort the data with index return
# s <- sort(Z, index.return = TRUE)
# L <- s$x
# sortix <- s$ix
# R <- I[sortix]
# 
# 
# # Turnbull estimator
# event <- !censored[sortix]
# # 3 for interval censoring and 0 for right censoring
# event[event==0 & is.finite(R)] <- 3
# type <- ifelse(all(is.finite(R)), "interval", "right")
# 
# if (type=="interval") {
#   fit  <- survfit(Surv(time=L, time2=R, event=event, type="interval") ~1, conf.type="plain", conf.int=1-alpha)
# } else {
#   fit  <- survfit(Surv(time=L, event=event, type="right") ~1, conf.type="plain", conf.int=1-alpha)
# }
# 
# f = stepfun(fit$time,c(1,fit$surv))
# est <- f(x)

