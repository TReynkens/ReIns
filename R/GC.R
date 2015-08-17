

# Gram-Charlier approximation for CDF
pGC <- function(x, moments = c(0,1,0,3), lower.tail = TRUE, log.p = FALSE) {
  
  L <- gcInput(moments)
  
  average <- L$average
  variance <- L$variance
  EZ3 <- L$EZ3
  EZ4 <- L$EZ4
  
  # Standardised x-values
  z <- (x-average)/sqrt(variance)
  
  # Gram-Charlier approximation
  p <- pnorm(z) + dnorm(z) * (-EZ3/6 * (z^2-1)  - (EZ4-3)/24 * (z^3-3*z))
  
  # Solve numerical issues to obtain valid probabilities
  p <- pmin(1, pmax(0,p))
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}


# Edgeworth approximation for CDF
pEdge <- function(x, moments = c(0,1,0,3), lower.tail = TRUE, log.p = FALSE) {
  
  L <- gcInput(moments)
  
  average <- L$average
  variance <- L$variance
  EZ3 <- L$EZ3
  EZ4 <- L$EZ4
  
  # Standardised x-values
  z <- (x-average)/sqrt(variance)
  
  # Gram-Charlier approximation
  p <- pnorm(z) + dnorm(z) * (-EZ3/6 * (z^2-1)  - (3*EZ4*(z^3-3*z) + EZ3^2*(z^5-10*z^3+15*z))/72 )
  
  # Solve numerical issues to obtain valid probabilities
  p <- pmin(1, pmax(0,p))
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}


# Check input and return (standardised) moments
gcInput <- function(moments) {
  
  if (length(moments)<4) {
    stop("Four moments should be provided.")
  }
  
  if (length(moments)>4) {
    warning("Only the first four moments are used.")
  }
  
  # Average
  average <- moments[1]
  # Variance
  variance <- moments[2] - moments[1]^2
  if (variance<=0) stop("Variance should be strictly positive.")
  
  # Standardised moment of order 3
  EZ3 <- standMoment(order=3, moments=moments)
  # Standardised moment of order 4
  EZ4 <- standMoment(order=4, moments=moments)
  
  if (EZ4<EZ3+1) {
    stop("The fourth standardised moments should be larger or equal to the 
         third standardised moment plus 1.")
  }

  return( list(average=average, variance=variance, EZ3=EZ3, EZ4=EZ4))
}

# Standardised (normalised) moment of order order
standMoment <- function(order = 1, moments) {
  
  m <- length(moments)
  if(order > m) {
    stop("The order should be smaller or equal to the length of moments.")
  }
  
  # First moment
  average <- moments[1]
  # Variance
  variance <- moments[2] - moments[1]^2
  
  standMom <- moments[order]
  
  # Add moment of order 0 (=1)
  moments <- c(1, moments)
  k <- 0:order
  # Subtract average and compute moment of order order
  # Use k+1 in last one because moments[1] is moment of order 0
  standMom <- sum(choose(order, k) * (-average)^(order-k) * moments[k+1])
  
  # Divide by standard deviation to the power order
  standMom <- standMom/sqrt(variance)^order
  
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

