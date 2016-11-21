
# This file contains the methods for approximating a CDF.

###########################################################################

# Classical approximations for the aggregate claim size distribution
pClas <- function(x, mean = 0, variance = 1, skewness = NULL, 
                  method = c("normal", "normal-power", "shifted Gamma", "shifted Gamma normal"), 
                  lower.tail = TRUE, log.p = FALSE) {
  
  
  # Check if valid input for method
  method <- match.arg(method)
  
  # Check input for mean and variance
  if (!is.numeric(mean) | !is.numeric(variance)) stop("Input arguments mean and variance should be numeric.")
  if (length(mean) > 1 | length(variance) > 1) stop("Input arguments mean and variance should have length 1.")
  if (variance <= 0) stop("The variance should be strictly positive.")
  
  # Check if skewness is provided
  if (is.null(skewness) & method != "normal") {
    stop(paste0("Input argument skewness cannot be NULL when using method \"", method, "\"."))
  }  
  
  # Check skewness
  if (!is.null(skewness)) {
    
    if (!is.numeric(skewness) | length(skewness) > 1) stop("skewness should be a numeric of length 1.")
    
    if (method != "normal" & skewness == 0) {
      warning(paste0("The skewness coefficient cannot be 0 when using \"", method, "\", 
                     the normal approximation will be used."))
      method <- "normal"
    }  
    
    if (method %in% c("shifted Gamma", "shifted Gamma normal") & skewness < 0) { 
      stop(paste0("The skewness coefficient should be strictly positive when using \"", method, "\"."))
    }

  }
  
  p <- numeric(length(x))
  
  if (method == "normal") {

    # Normal approximation using mean and variance
    p <- pnorm((x - mean) / sqrt(variance))
    
  } else if (method == "normal-power") {
    
    # Normal-power approximation: correction of normal approximation using skewness coefficient
    
    z <- (x - mean) / sqrt(variance)
    # Problems with z<1
    if (any(z < 1)) {
      ind <- which(z >= 1)
      warning("Only estimates for F(x) for values of x larger than or equal to mean + sqrt(variance) are provided.")
      p[-ind] <- NaN
      
    } else {
      ind <- 1:length(x)
    }

    # Problems with negative values in root
    if (any( (9 / skewness ^ 2 + 6 * z[ind] / skewness + 1) < 0)) {
      ind2 <- ind[which( (9 / skewness ^ 2 + 6 * z[ind] / skewness + 1) >= 0)]
      warning("Only estimates for F(x) for values of x where \'9/nu^2 + 6*(x-mu)/sigma/nu + 1 > 0\' are provided.")
      p[-ind2] <- NaN
      
    } else {
      ind2 <- ind
    }
    
    p[ind2] <- pnorm(sqrt(9 / skewness ^ 2 + 6 * z[ind2] / skewness + 1) - 3 / skewness)
    
  } else if (method == "shifted Gamma") {
    
    # Shifted Gamma approximation

    y <- x - mean + 2 * sqrt(variance) / skewness
    
    p <- pgamma(y, shape = 4 / skewness ^ 2, rate = 2 / (skewness * sqrt(variance)))
    
  } else {
    
    # Normal approximation to shifted Gamma distribution
    
    z <- (x - mean) / sqrt(variance)
    # Problems with z<1
    if (any( z < 1 )) {
      ind <- which(z >= 1)
      warning("Only estimates for F(x) for values of x larger than or equal to mean + sqrt(variance) are provided.")
      p[-ind] <- NaN
      
    } else {
      ind <- 1:length(x)
    }
    
    p[ind] <- pnorm(sqrt(16 / skewness ^ 2 + 8 * z[ind] / skewness) - sqrt(16 / skewness ^ 2 - 1))
  }
  
  # Solve numerical issues to obtain valid probabilities
  p <- pmin(1, pmax(0, p))
  
  if (!lower.tail) p <- 1 - p
  
  if (log.p) p <- log(p)
  
  return(p)
  
}


# Gram-Charlier approximation for CDF
pGC <- function(x, moments = c(0, 1, 0, 3), raw = TRUE, lower.tail = TRUE, log.p = FALSE) {
  
  if (length(moments) < 4) {
    stop("Four moments should be provided.")
  }
  
  if (length(moments) > 4) {
    warning("Only the first four moments are used.")
  }
  
  # Use raw moments of X
  if (raw) {

    # Mean
    average <- moments[1]
    # Variance
    variance <- moments[2] - moments[1] ^ 2
    
    # Standardised moment of order 3
    EZ3 <- .standMoment(order = 3, moments = moments)
    # Standardised moment of order 4
    EZ4 <- .standMoment(order = 4, moments = moments)
 
  } else {
    # Mean, variance, skewness and kurtosis are provided
    average <- moments[1]
    variance <- moments[2]
    EZ3 <- moments[3]
    EZ4 <- moments[4]
    
  }
  
  if (variance <= 0) stop("The variance should be strictly positive.")
  
  if (EZ4 < EZ3 + 1) {
    stop("The fourth standardised moment should be larger than or equal to the 
         third standardised moment plus 1.")
  }

  # Standardised x-values
  z <- (x - average) / sqrt(variance)
  
  # Gram-Charlier approximation
  p <- pnorm(z) + dnorm(z) * (- EZ3 / 6 * (z ^ 2 - 1)  - (EZ4 - 3) / 24 * (z ^ 3 - 3 * z))
  
  # Solve numerical issues to obtain valid probabilities
  p <- pmin(1, pmax(0, p))
  
  if (!lower.tail) p <- 1 - p
  
  if (log.p) p <- log(p)
  
  return(p)
}


# Edgeworth approximation for CDF
pEdge <- function(x, moments = c(0, 1, 0, 3), raw = TRUE, lower.tail = TRUE, log.p = FALSE) {
  
  
  if (length(moments) < 4) {
    stop("Four moments should be provided.")
  }
  
  if (length(moments) > 4) {
    warning("Only the first four moments are used.")
  }
  
  # Use raw moments of X
  if (raw) {

    # Mean
    average <- moments[1]
    # Variance
    variance <- moments[2] - moments[1] ^ 2
    
    # Standardised moment of order 3
    EZ3 <- .standMoment(order = 3, moments = moments)
    # Standardised moment of order 4
    EZ4 <- .standMoment(order = 4, moments = moments)
    
  } else {
    # Mean, variance, skewness and kurtosis are provided
    average <- moments[1]
    variance <- moments[2]
    EZ3 <- moments[3]
    EZ4 <- moments[4]
    
  }
  
  if (variance <= 0) stop("The variance should be strictly positive.")
  
  if (EZ4 < EZ3 + 1) {
    stop("The fourth standardised moment should be larger than or equal to the 
         third standardised moment plus 1.")
  }
  
  # Standardised x-values
  z <- (x - average) / sqrt(variance)
  
  # Edgeworth approximation
  p <- pnorm(z) + dnorm(z) * (- EZ3 / 6 * (z ^ 2 - 1)  - (3 * EZ4 * (z ^ 3 - 3 * z) + 
                               EZ3 ^ 2 * (z ^ 5 - 10 * z ^ 3 + 15 * z)) / 72 )
  
  # Solve numerical issues to obtain valid probabilities
  p <- pmin(1, pmax(0, p))
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}



# Standardised (normalised) moment of order order
.standMoment <- function(order = 1, moments) {
  
  m <- length(moments)
  if (order > m) {
    stop("The order should be smaller than or equal to the length of moments.")
  }
  
  # First moment
  average <- moments[1]
  # Variance
  variance <- moments[2] - moments[1] ^ 2
  
  standMom <- moments[order]
  
  # Add moment of order 0 (=1)
  moments <- c(1, moments)
  k <- 0:order
  # Subtract average and compute moment of order order
  # Use k+1 in last one because moments[1] is moment of order 0
  standMom <- sum(choose(order, k) * (- average) ^ (order - k) * moments[k + 1])
  
  # Divide by standard deviation to the power order
  standMom <- standMom / sqrt(variance) ^ order
  
  return(standMom)
}



# # Gram-Charlier approximation for density
# dGC <- function(x, moments = c(0,1,0,3), log = FALSE) {
#   
#   L <- gcInput(moments)
#   
#   average <- L$average
#   variance <- L$variance
#   EZ3 <- L$EZ3
#   EZ4 <- L$EZ4
#   
#   # Standardised x-values
#   z <- (x-average)/sqrt(variance)
#   
#   # Gram-Charlier approximation
#   d <- dnorm(z) * (1 + EZ3/6 * (z^3-3*z) + (EZ4-3)/24 * (z^4-6*z^2+3))
#   
#   # Solve numerical issues
#   d <- pmax(0,d)
#   
#   
#   if (log) d <- log(d)
#   
#   return(d)
# }
# 
# # Edgeworth approximation for density
# dEdge <- function(x, moments = c(0,1,0,3), log = FALSE) {
#   
#   L <- gcInput(moments)
#   
#   average <- L$average
#   variance <- L$variance
#   EZ3 <- L$EZ3
#   EZ4 <- L$EZ4
#   
#   # Standardisesd x-values
#   z <- (x-average)/sqrt(variance)
#   
#   # Gram-Charlier approximation
#   d <- dnorm(z) * ( 1 + EZ3/6 * (z^2-2*z-1)  + (3*EZ4*(z^4-6*z^2+3) + EZ3^2*(z^6+15*z^4+25*z^2-15))/72 )
#   
#   # Solve numerical issues to obtain valid probabilities
#   d <- pmax(0,d)
#   
#   if (log) d <- log(d)
#   
#   return(d)
# }

