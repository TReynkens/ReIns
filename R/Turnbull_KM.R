
# This file contains the functions for the Kaplan-Meier and Turnbull estimators.

###########################################################################


# Kaplan-Meier estimator for the CDF evaluated in x
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary sample CDF
KaplanMeier <- function(x, data, censored, conf.type = "plain", conf.int = 0.95) {

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
  
  fit  <- survfit(Surv(X, delta) ~1, conf.type=conf.type, conf.int=conf.int)

  f <- stepfun(fit$time, c(1, fit$surv))
  est <- f(x)
  
  return(list(surv = est, fit=fit))
}


# Turnbull estimatorfor the survival function evaluated in x
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
  if (length(L) != length(R)) {
    stop("L and R should have the same length.")
  }
  
  if (length(trunclower) != 1) {
    stop("trunclower should have length 1.")
  }
  
  if (length(truncupper) != 1) {
    stop("truncupper should have length 1.")
  }
  
  if (any(L < trunclower)) {
    stop("All elements of L should be larger than trunclower.")
  }
  
  if (any(R > truncupper)) {
    stop("All elements of R should be smaller than truncupper.")
  }
  
  if (any(L > R)) {
    stop("Each element of L should be smaller or equal than each corresponding element of R.")
  }
  
  if (!is.numeric(x)) {
    stop("x should be numeric.")
  }
  
  censored <- .checkCensored(censored, length(L))
  
  
  ft <- .Turnbull_internal(L=L, R=R, censored=censored, conf.type = conf.type, conf.int = conf.int,
                           trunclower=trunclower, truncupper=truncupper)
  est <- ft$f(x)
  
  # Turnbull estimator for the survival function
  return( list(surv = est, fit=ft$fit) )
}


# Function to return Turnbull stepfunction for survival function
.Turnbull_internal <- function(L, R, censored, trunclower = 0, truncupper = Inf, 
                               conf.type = "plain", conf.int = 0.95) {
  
  # Uncensored observations
  R[!censored] <- L[!censored]
  # Right censored observations, set upper bound to NA
  R[censored & R == truncupper] <- NA
  # Left censored observations, set lower bound to NA
  L[censored & L == trunclower] <- NA


  fit  <- survfit(Surv(time=L, time2=R, type="interval2") ~1,
                  conf.type = conf.type, conf.int = conf.int)

  # Survival function (right continuous)
  f <- stepfun(fit$time, c(1, fit$surv), right=FALSE)
  
  # Quantile function (left continuous)
  fquant <- function(p) {
    
    fq <- stepfun(1-fit$surv, c(fit$time, Inf), right=TRUE)
    
    return(ifelse(p >= 0 & p <= 1, fq(p), NaN))
  }

  return(list(f=f, fit=fit, fquant=fquant))
}


# Turnbull estimator for the survival function using interval package evaluated in x
# L and R are the lower and upper values for interval censoring
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary sample CDF
#
.Turnbull2 <- function(x, L, R, censored, trunclower = 0, truncupper = Inf) {
  
  
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
  if (length(L) != length(R)) {
    stop("L and R should have the same length.")
  }
  
  if (length(trunclower) != 1) {
    stop("trunclower should have length 1.")
  }
  
  if (length(truncupper) != 1) {
    stop("truncupper should have length 1.")
  }
  
  if (any(L < trunclower)) {
    stop("All elements of L should be larger than trunclower.")
  }
  
  if (any(R > truncupper)) {
    stop("All elements of R should be smaller than truncupper.")
  }
  
  if (any(L > R)) {
    stop("Each element of L should be smaller or equal than each corresponding element of R.")
  }
  
  if (!is.numeric(x)) {
    stop("x should be numeric.")
  }
  censored <- .checkCensored(censored, length(L))
  
  
  ft <- .Turnbull_internal2(L=L, R=R, censored=censored, trunclower=trunclower, truncupper=truncupper)
  est <- ft$f(x)
  
  # Turnbull estimator for the survival function
  return( list(surv = est, fit=ft$fit) )
}

# Function to return Turnbull step function for survival function using interval package
.Turnbull_internal2 <- function(L, R, censored, trunclower = 0, truncupper = Inf) {
  
  # Uncensored observations
  R[!censored] <- L[!censored]
  # Right censored observations, set upper bound to Inf
  R[censored & R == truncupper] <- Inf
  # Left censored observations, set lower bound to 0
  L[censored & L == trunclower] <- 0
  
  # Fit Turnbull estimator using interval package,
  # starting value using Icens package, use no starting value if problems
  fit <- tryCatch(interval::icfit(Surv(time=L, time2=R, type="interval2")~1, 
                                  conf.int=FALSE, initfit=Icens::EMICM(cbind(L, R))),
                  error = function(e) interval::icfit(Surv(time=L, time2=R, type="interval2")~1, 
                                                      conf.int=FALSE, initfit=NULL))
  
  # Extract Turnbull intervals
  x <- fit$intmap
  # Extract probabilities corresponding to these intervals
  y <- fit$pf
  # Only keep intervals where probabilities are >0
  x <- x[,y>0]
  y <- y[y>0]

  # Make function with linear interpolation in Turnbull intervals
  # Plot empirical survival function with linear interpolation in TB intervals
  xl <- x[1,]
  xu <- x[2,]
  
  # Left and right probabilities
  probl <- c(1, 1-cumsum(y)[-length(y)])
  probu <- 1-cumsum(y)
  
  # Correct names of probl
  names(probl) <- names(probu)
  
  # Subtract small number to make function right continuous
  #xl <- xl - sqrt(.Machine$double.eps)
  eps <- sqrt(.Machine$double.eps)
  ind_jump <- which(abs(xu-xl)<eps)
  xl[ind_jump] <- xl[ind_jump] - eps
  
  survall <- c(probl, probu)
  xall <- c(xl, xu)
  
  # Order x and y-values
  survall <- survall[order(xall)]
  xall <- xall[order(xall)]

  # Linear interpolation with jumps in ties
  f <- approxfun(xall, survall, ties = "ordered", rule=2, yleft=1)

  ##
  # Make quantile function
  
  xl <- x[1,]
  xu <- x[2,]
  
  # Left and right probabilities
  probl <- c(0, cumsum(y)[-length(y)])
  probu <- cumsum(y)
  
  # Correct names of probl
  names(probl) <- names(probu)
  
  # Add small number to make function left continuous
  probl <- probl + 10^(-10)
  
  # Combine values of x
  xall2 <- c(xl, xu)
  
  # Combine probabilities
  proball <- c(probl, probu)
  
  # Order x and y-values
  xall2 <- xall2[order(proball)]
  proball <- proball[order(proball)]
  
  # Quantile function (left continuous)
  fquant <- function(p) {
    
    fq <- approxfun(proball, xall2, ties = "ordered", rule=2)
    
    return(ifelse(p >= 0 & p <= 1, fq(p), NaN))
  }
  
  return(list(f=f, fit=fit, xall=xall, survall=survall, fquant=fquant))
}

