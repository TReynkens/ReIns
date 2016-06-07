

# Mean excess plot using Turnbull estimator
MeanExcess_TB <- function(L, U = L, censored, trunclower = 0, truncupper = Inf, 
                          plot = TRUE, k = TRUE, main = "Mean excess plot", ...) {
  
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
  me <- numeric(n)
  me2 <- numeric(n)
  
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
    
    if (k) {
      .plotfun(K, me[K], type="p", xlab="k", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    } else {
      .plotfun(ic[n-K], me[K], xlab=bquote(X["n-k,n"]), type="p", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    }
    
  }
  
  # Output list with values of k, order statistics X_n-k,n 
  # and mean excess scores e_k,n
  .output(list(k=K, X=ic[n-K], e=me[K]), plot=plot, add=FALSE)
}


