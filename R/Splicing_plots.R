
# This file contains all plots to assess the goodness-of-fit of the splicing model.

###########################################################################


# Plot of fitted survival function and ECDF estimator + bounds
SpliceECDF <- function(x, X, splicefit, alpha = 0.05, ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  plot(x, 1-pSplice(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)
  
  # ECDF estimator
  fit  <- ecdf(X)
  # Empirical surivival function
  est <- 1-fit(x)
  
  # Confidence bounds using Dvoretzky-Kiefer-Wolfowitz inequality
  # http://stats.stackexchange.com/questions/55500/confidence-intervals-for-empirical-cdf
  n <- length(X)
  eps <- sqrt(1/(2*n)*log(2/alpha))
  lines(x, est, lty=1, col="red")
  lines(x, pmax(est-eps, 0), lty=2, col="blue")
  lines(x, pmin(est+eps, 1), lty=2, col="blue")
  legend("topright", c("Fitted survival function", "Empirical survival function", "95% confidence bounds"),
         lty=c(1, 1, 2), col=c("black", "red", "blue"))
}

# Plot of fitted survival function and Turnbull estimator + confidence intervals
SpliceTB <- function(x = sort(L), L, U = L, censored, splicefit, alpha = 0.05, ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # Check lengths
  if (length(L) != length(U)) stop("L and U should have the same length.")
  
  if (length(censored) == 1) censored <- rep(censored, length(L))
  
  if (length(censored) != length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Plot fitted survival function
  plot(x, 1-pSplice(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)
  
  
  # Add Turnbull survival function
  tb <- Turnbull(x, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                 truncupper=max(splicefit$EVTfit$endpoint), conf.type="plain", conf.int=1-alpha)
  lines(x, tb$surv, col="red")
  
  # Add confidence intervals
  fit <- tb$fit
  lines(fit$time, fit$lower, col = "blue", lty = 2)
  lines(fit$time, fit$upper, col = "blue", lty = 2)
  legend("topright", c("Fitted survival function", "Turnbull estimator", "95% confidence intervals"),
         lty=c(1, 1, 2), col=c("black", "red", "blue"))
}


# Probability - probability plot with ECDF
SplicePP <- function(X, splicefit, x = sort(X), log = FALSE, plot = TRUE, main = "Splicing PP-plot", ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # ECDF estimator
  fit  <- ecdf(X)
  est <- 1-fit(x)
  
  # Plot fitted survival function vs. empirical survival function or use minus log-versions
  if (log) {
    ind <- est>0
    spp.the <- -log(1-pSplice(x[ind], splicefit=splicefit))
    spp.emp <- -log(est[ind])
    .plotfun(spp.emp, spp.the, type="p",
             xlab="-log(Empirical survival probability)",
             ylab="-log(Fitted survival probability)", plot=plot, add=FALSE, main=main, ...)
    
  } else {
    spp.the <- 1-pSplice(x, splicefit=splicefit)
    spp.emp <- est
    .plotfun(spp.emp, spp.the, type="p", xlab="Empirical survival probability",
             ylab="Fitted survival probability", plot=plot, add=FALSE, main=main, ...)
  }
  # Add 45 degree line
  if (plot) abline(0, 1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(spp.the=spp.the, spp.emp=spp.emp), plot=plot, add=FALSE)
}

# Probability - probability plot with Turnbull estimator
SplicePP_TB <- function(L, U = L, censored, splicefit, x = NULL, log = FALSE, plot = TRUE, main = "Splicing PP-plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.null(x) & !is.numeric(x)) stop("x should be a numeric vector or NULL.")
  
  # Check lengths
  if (length(L) != length(U)) stop("L and U should have the same length.")
  
  if (length(censored) == 1) censored <- rep(censored, length(L))
  
  if (length(censored) != length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                  truncupper=max(splicefit$EVTfit$endpoint))
    if (is.null(x)) {
      # # Use unique points of Turnbull intervals since Turnbull estimator is exact there
      # x <- unique(as.numeric(SurvTB$fit$intmap))
      # Use empirical quantiles
      n <- length(L)
      p <- (1:n) / (n+1)
      x <- SurvTB$fquant(p)
    }
    
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                 truncupper=max(splicefit$EVTfit$endpoint))
    
    if (is.null(x)) {
      # # Use knots (jump points)
      # x <- unique(SurvTB$fit$knots)
      # Use empirical quantiles
      n <- length(L)
      p <- (1:n) / (n+1)
      x <- SurvTB$fquant(p)
    }
  }
  
  # Survival function in x
  surv <- SurvTB$f(x)
  
  # Plot fitted survival function vs. Turnbull survival function or use minus log-versions
  if (log) {
    ind <- surv>0
    spp.the <- -log(1-pSplice(x[ind], splicefit=splicefit))
    spp.emp <- -log(surv[ind])
    .plotfun(spp.emp, spp.the, type="p",
             xlab="-log(Turnbull survival probability)",
             ylab="-log(Fitted survival probability)", plot=plot, add=FALSE, main=main, ...)
    
  } else {
    spp.the <- 1-pSplice(x, splicefit=splicefit)
    spp.emp <- surv
    .plotfun(spp.emp, spp.the, type="p", xlab="Turnbull survival probability",
             ylab="Fitted survival probability", plot=plot, add=FALSE, main=main, ...)
  }
  # Add 45 degree line
  if (plot) abline(0, 1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(spp.the=spp.the, spp.emp=spp.emp), plot=plot, add=FALSE)
}



# Quantile - quantile plot with empirical quantiles
SpliceQQ <- function(X, splicefit, p = NULL, plot = TRUE, main = "Splicing QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(X, pos=FALSE)
  
  X <- as.numeric(sort(X))
  n <- length(X)
  
  # Probabilities
  if (is.null(p)) {
    p <- (1:n)/(n+1)
    
  } else {
    
    if (!is.numeric(p)) {
      stop("p should be numeric.")
    }
    
    if (is.infinite(max(splicefit$EVTfit$endpoint)) & any(p == 1)) {
      stop("All elements of p should be strictly smaller than 1 since the splicing
           distribution has an infinite endpoint.")
    }
    
    if (any(p < 0 | p > 1)) {
      stop("All elements of p should be between 0 and 1.")
    }
    
    if (length(p) != n) {
      stop("p should have the same length as x.")
    }
    
    p <- sort(p)
    }
  
  
  # Empirical quantiles
  sqq.emp <- sort(X)
  
  # Quantiles of fitted distribution
  sqq.the <- qSplice(p=p, splicefit=splicefit)
  
  .plotfun(sqq.the, sqq.emp, type="p", xlab="Quantiles of splicing fit", ylab="Empirical quantiles", 
           main=main, plot=plot, add=FALSE, ...)
  
  # Add 45 degree line
  if (plot) abline(0, 1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(sqq.the=sqq.the, sqq.emp=sqq.emp), plot=plot, add=FALSE)
}


# QQ-plot with Turnbull estimator
SpliceQQ_TB <- function(L, U = L, censored, splicefit, p = NULL, plot = TRUE, main = "Splicing QQ-plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  
  # Check lengths
  if (length(L) != length(U)) stop("L and U should have the same length.")
  
  if (length(censored) == 1) censored <- rep(censored, length(L))
  
  if (length(censored) != length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE) & !all(censored == rep(0, length(L)))) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                  truncupper=max(splicefit$EVTfit$endpoint))
    
    # # Use unique points of Turnbull intervals since Turnbull estimator is exact there
    # x <- unique(as.numeric(SurvTB$fit$intmap))
    # # Corresponding function values
    # s <- SurvTB$f(x)
    # p <- 1-s
    
  } else {
    
    # Special warning when no censoring
    if (all(censored == rep(0, length(L)))) {
      warning("Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    } else {
      warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    }
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                 truncupper=max(splicefit$EVTfit$endpoint))
    
    # # Use knots (jump points)
    # x <- knots(SurvTB$f)
    # # Corresponding function values
    # s <- eval(expression(y), envir = environment(SurvTB$f))
    # p <- 1-s
  }
  
  
  # Probabilities
  n <- length(L)
  if (is.null(p)) {
    p <- (1:n)/(n+1)
    
  } else {
    
    if (!is.numeric(p)) {
      stop("p should be numeric.")
    }
    
    if (is.infinite(max(splicefit$EVTfit$endpoint)) & any(p == 1)) {
      stop("All elements of p should be strictly smaller than 1 since the splicing
           distribution has an infinite endpoint.")
    }
    
    if (any(p < 0 | p > 1)) {
      stop("All elements of p should be between 0 and 1.")
    }
    
    if (length(p) != n) {
      stop("p should have the same length as x.")
    }
    
    p <- sort(p)
    }
  # Use empirical quantiles
  x <- SurvTB$fquant(p)
  
  # Remove 1 if infinite endpoint
  if (is.infinite(max(splicefit$EVTfit$endpoint))) {
    sqq.emp <- x[p<1-10^(-5)]
    p <- p[p<1-10^(-5)]
  } else {
    sqq.emp <- x
  }
  
  # Quantiles of fitted distribution
  sqq.the <- qSplice(p=p, splicefit=splicefit)
  
  .plotfun(sqq.the, sqq.emp, type="p", xlab="Quantiles of splicing fit", ylab="Turnbull quantiles", 
           main=main, plot=plot, add=FALSE, ...)
  
  # Add 45 degree line
  if (plot) abline(0, 1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(sqq.the=sqq.the, sqq.emp=sqq.emp), plot=plot, add=FALSE)
}





# Log-log plot with empirical survival function and fitted survival function
SpliceLL <- function(x = sort(X), X, splicefit, plot = TRUE, main = "Splicing LL-plot", ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # ECDF estimator
  fit  <- ecdf(X)
  X <- sort(X)
  
  sll.emp <- log(1-fit(X))
  sll.the <- log(1-pSplice(x, splicefit=splicefit))
  
  if (plot) {
    # Plot log of empirical survival function vs. sorted values
    plot(log(X), sll.emp, ylab="log(empirical survival probability)", xlab="log(X)", type="p", main=main, ...)
    # Add log of fitted survival function
    lines(log(x), sll.the)
  }
  
  .output(list(logX=log(X), sll.emp=sll.emp, logx=log(x), sll.the=sll.the), plot=plot, add=FALSE)
}


# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB <- function(x = sort(L), L, U = L, censored, splicefit, plot = TRUE, main = "Splicing LL-plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # Check lengths
  if (length(L) != length(U)) stop("L and U should have the same length.")
  
  if (length(censored) == 1) censored <- rep(censored, length(L))
  
  if (length(censored) != length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  
  Zs <- sort(L)
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull2(Zs, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                         truncupper=max(splicefit$EVTfit$endpoint))
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    SurvTB <- Turnbull(Zs, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                       truncupper=max(splicefit$EVTfit$endpoint))
  }
  
  sll.emp <- log(SurvTB$surv)
  sll.the <- log(1-pSplice(x, splicefit=splicefit))
  
  if (plot) {
    # Plot log of Turnbull survival function vs. sorted values
    plot(log(Zs), sll.emp, ylab="log(Turnbull survival probability)", xlab="log(X)", type="p", main=main, ...)
    # Add log of fitted survival function
    lines(log(x), sll.the)
  }
  .output(list(logX=log(Zs), sll.emp=sll.emp, logx=log(x), sll.the=sll.the), plot=plot, add=FALSE)
}
