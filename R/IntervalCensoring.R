

# Mean excess plot using Turnbull estimator
MeanExcess_TB <- function(L, U = L, censored, trunclower = 0, truncupper = Inf, 
                          plot = TRUE, k = FALSE, intervalpkg = TRUE, main = "Mean excess plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  kplot <- k
  
  # Turnbull survival function
  # Use interval package if available and if requested by user
  if (requireNamespace("interval", quietly = TRUE) & intervalpkg) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$xall
    y <- SurvTB$survall
    m <- length(x)
    
  } else {
    # Issue warning if interval package requested by user but not available
    if (intervalpkg) {
      warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    }

    
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$fit$time
    y <- SurvTB$fit$surv
    m <- length(x)
  }
  
  
  # Values to evaluate mean excess function in
  n <- length(L)
  p <- (1:n) / (n+1)
  ic <- SurvTB$fquant(p)
  
  # Keep unique elements
  eps <- sqrt(.Machine$double.eps)
  ic[which(diff(ic)<=eps)+1] <- NA
  ic <- ic[!is.na(ic)]
  
  # Remove infinite elements
  ic <- ic[is.finite(ic)]
  
  
  n <- length(ic)
  K <- 1:(n-1)
  me <- numeric(n)

  # Compute slopes and intercepts of function in intervals
  as <- diff(y) / diff(x)
  as[diff(x)<=eps] <- 0
  bs <- y[1:(m-1)] - as * x[1:(m-1)]
  # Remove influence of jumps
  bs[is.infinite(as)] <- 0
  as[is.infinite(as)] <- 0
  
  for (k in K) {
    
    # Threshold
    r <- ic[n-k]
    
    if (r>x[m]) {
      me[k] <- 0
      
    } else {
      
      # Find correct interval (right-closed)
      ind_r <- length(x) - findInterval(-r, -rev(x))
      #ind_r = findInterval(r, x)
      
      if (r<x[1]) r <- x[1]
      
      # Slope and intercept of function in interval containing r
      if(ind_r==0) {
        # Special case if ind_r <- 0
        a <- 0
        b <- 1
        
      } else {
        a <- as[ind_r]
        b <- bs[ind_r]
      }
      
      # Use the fact that the Turnbull survival function is piecewise linear with jumps
      # First, contribution of interval containing r, then contributions of all intervals after this interval
      me[k] <- a * (x[ind_r+1]^2-r^2) / 2 + b *(x[ind_r+1]-r)  + 
        sum(as[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]^2) / 2 + bs[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]))
    }
    
  }
  
  # Mean-excess is int_r^Inf (1-F(z)) dz / (1-F(r))
  me[K] <- me[K] / SurvTB$f(ic[n-K])
  
  # Plot estimates
  if (plot) {
    
    if (kplot) {
      .plotfun(K, me[K], type="p", xlab="k", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    } else {
      .plotfun(ic[n-K], me[K], xlab=bquote(X["n-k,n"]), type="p", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    }
    
  }
  
  # Output list with values of k, order statistics X_n-k,n 
  # and mean excess scores e_k,n
  .output(list(k=K, X=ic[n-K], e=me[K]), plot=plot, add=FALSE)
}


# Pareto QQ-plot adapted for interval censoring using Turnbull estimator
icParetoQQ <- function(L, U = L, censored, trunclower=0, truncupper=Inf, plot = TRUE, main = "Pareto QQ-plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE) & !all(censored==rep(0, length(L)))) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=trunclower,
                                  truncupper=truncupper)
    
  } else {
    
    # Special warning when no censoring
    if (all(censored==rep(0, length(L)))) {
      warning("Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    } else {
      warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    }
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=trunclower,
                                 truncupper=truncupper)
  }
  
  
  # Probabilities
  n <- length(L)
  p <- (1:n) / (n+1)
  # Empirical quantiles
  x <- SurvTB$fquant(p)
  
  # Remove 1 if infinite endpoint
  if(is.infinite(truncupper)) {
    pqq.emp <- log(x)[p<1-10^(-5)]
    p <- p[p<1-10^(-5)]
  } else {
    pqq.emp <- log(x)
  }
  
  
  # Quantiles of fitted distribution
  pqq.the <- -log(1 - p)
  
  .plotfun(pqq.the, pqq.emp, type="p", xlab="Quantiles of standard exponential", ylab="Log. of empirical quantiles", 
           main=main, plot=plot, add=FALSE, ...)
  
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(pqq.the=pqq.the, pqq.emp=pqq.emp), plot=plot, add=FALSE)
}

################################################################################################

# Estimator for gamma using Turnbull estimator
icHill <- function(L, U, censored, trunclower = 0, truncupper = Inf, 
                   logk = FALSE, plot = TRUE, add = FALSE, main = "Hill estimates of EVI", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$xall
    y <- SurvTB$survall
    m <- length(x)
    
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$fit$time
    y <- SurvTB$fit$surv
    m <- length(x)
  }

  # Values to evaluate mean excess function in
  n <- length(L)
  p <- (1:n) / (n+1)
  ic <- SurvTB$fquant(p)
  
  # Keep unique elements
  eps <- sqrt(.Machine$double.eps)
  ic[which(diff(ic)<=eps)+1] <- NA
  ic <- ic[!is.na(ic)]
  
  
  n <- length(ic)
  K <- 1:(n-1)
  gamma <- numeric(n)
  
  
  # Compute slopes and intercepts of function in intervals
  as <- diff(y) / diff(x)
  as[diff(x)<=eps] <- 0
  bs <- y[1:(m-1)] - as * x[1:(m-1)]
  # Remove influence of jumps
  bs[is.infinite(as)] <- 0
  as[is.infinite(as)] <- 0
  
  
  for (k in K) {
    
    # Threshold
    r <- ic[n-k]
    
    if (r>x[m]) {
      gamma[k] <- 0
      
    } else {
      
      # Find correct interval (right-closed)
      ind_r <- length(x) - findInterval(-r, -rev(x))
      #ind_r = findInterval(r, x)
      
      if (r<x[1]) r <- x[1]
      
      
      # Slope and intercept of function in interval containing r
      if(ind_r==0) {
        # Special case if ind_r <- 0
        a <- 0
        b <- 1
        
      } else {
        a <- as[ind_r]
        b <- bs[ind_r]
      }
      
      # Use the fact that the Turnbull survival function is piecewise linear with jumps
      # First, contribution of interval containing r, then contributions of all intervals after this interval
      gamma[k] <- a * (x[ind_r+1]-r) + b *(log(x[ind_r+1])-log(r))  + 
        sum(as[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]) + bs[(ind_r+1):(m-1)] * diff(log(x[(ind_r+1):m])))
    }
  }
  
  # Estimate for gamma is int_r^Inf (1-F(z))/z dz / (1-F(r))
  gamma[K] <- gamma[K] / SurvTB$f(ic[n-K])
  
  # Plot estimates
  if (logk) {
    .plotfun(log(K), gamma[K], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, gamma[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }

  # Output list with values of k, estimates for gamma
  # and order statistics X_n-k,n 
  .output(list(k=K, gamma=gamma[K], X=ic[n-K]), plot=plot, add=FALSE)
}


# Estimate for gamma using Turnbull estimator using a given threshold
.icHill_t <- function(L, U, censored, threshold, trunclower = 0, truncupper = Inf) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$xall
    y <- SurvTB$survall
    m <- length(x)
    
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$fit$time
    y <- SurvTB$fit$surv
    m <- length(x)
  }
  
  # Values to evaluate mean excess function in
  n <- length(L)
  p <- (1:n) / (n+1)
  ic <- SurvTB$fquant(p)
  
  # Keep unique elements
  eps <- sqrt(.Machine$double.eps)
  ic[which(diff(ic)<=eps)+1] <- NA
  ic <- ic[!is.na(ic)]

  # Compute slopes and intercepts of function in intervals
  as <- diff(y) / diff(x)
  as[diff(x)<=eps] <- 0
  bs <- y[1:(m-1)] - as * x[1:(m-1)]
  # Remove influence of jumps
  bs[is.infinite(as)] <- 0
  as[is.infinite(as)] <- 0
  
  # Threshold
  r <- threshold
  
  if (r>x[m]) {
    gamma <- 0
    
  } else {
    
    # Find correct interval (right-closed)
    ind_r <- length(x) - findInterval(-r, -rev(x))
    #ind_r = findInterval(r, x)
    
    if (r<x[1]) r <- x[1]
    
    
    # Slope and intercept of function in interval containing r
    if(ind_r==0) {
      # Special case if ind_r <- 0
      a <- 0
      b <- 1
      
    } else {
      a <- as[ind_r]
      b <- bs[ind_r]
    }
    
    # Use the fact that the Turnbull survival function is piecewise linear with jumps
    # First, contribution of interval containing r, then contributions of all intervals after this interval
    gamma <- a * (x[ind_r+1]-r) + b *(log(x[ind_r+1])-log(r))  + 
      sum(as[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]) + bs[(ind_r+1):(m-1)] * diff(log(x[(ind_r+1):m])))
  }
 
  # Estimate for gamma is int_r^Inf (1-F(z))/z dz / (1-F(r))
  gamma <- gamma / SurvTB$f(threshold)
  
  return(as.numeric(gamma))
}
