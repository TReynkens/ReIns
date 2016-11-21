
# This file contains implementations of the risk measures: excess-loss premiums (for the splicing model
# and several EVT models), and Value-at-Risk (VaR) and Expected Shortfall (ES) for the splicing model.

###########################################################################

# Integrated tail function of (truncated) Pareto distribution using Hill estimates
# Pi(R)
.IntTailPareto <- function(data, gamma, R, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                    main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  .checkInput(data, gamma)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  if (length(R)>1 & length(R) != (n-1)) {
    stop("R should be a numeric of length 1 or n-1.")
  }
  
  # Premium
  premium[K] <- (K+1) / (n+1) * .IntTailPareto_aux(t=X[n-K], gamma=gamma[K], R=R, endpoint=endpoint)
  
  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when R >= X[n-K]
  premium[X[n-K]>R] <- NA
  
  if (any(X[n-K]>R) & warnings) {
    warning("R is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R), plot=plot, add=add)
  
}

# Auxiliary function used in ExcessSplice
.IntTailPareto_aux <- function(t, gamma, R, endpoint = Inf) {
  
  if (any(gamma >= 1)) stop("gamma should be strictly smaller than 1.")
  
  # Premium
  if (is.finite(endpoint)) {

    # Truncated Pareto model
    beta <- endpoint/t
    premium <- 1 / (1-beta^(-1/gamma)) * ( t^(1/gamma) / (1/gamma-1) * 
                                             (R^(1-1/gamma) - endpoint^(1-1/gamma)) + beta^(-1/gamma) * (R - endpoint) )
    # Special case if R>=endpoint
    premium[R >= endpoint] <- 0
    
  } else {
    # Pareto model
    premium <- 1 / (1/gamma-1) * t^(1/gamma) * R^(1-1/gamma)
  }
  return(premium)
}


# Integrated tail function using EPD estimates
# Pi(R)
.IntTailEPD <- function(data, gamma, kappa, tau, R, warnings = TRUE, plot = TRUE, add = FALSE,
                    main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  if (length(R) > 1) {
    stop("R should be a numeric of length 1.")
  }
  
  if (length(gamma) != length(kappa) | length(gamma) != length(tau)) {
    stop("gamma, kappa and tau should have equal length.")
  } 
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  premium[K] <- (K+1) / (n+1) * X[n-K]^(1/gamma[K]) * ( (1-kappa[K]/gamma[K]) / (1/gamma[K]-1) * R^(1-1/gamma[K]) +
                                                          kappa[K]/(gamma[K]*X[n-K]^tau[K]) / (1/gamma[K]-tau[K]-1) *
                                                          R^(1+tau[K]-1/gamma[K]) )
  
  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when R >= X[n-K]
  premium[X[n-K]>R] <- NA
  
  if (any(X[n-K]>R) & warnings) {
    warning("R is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R), plot=plot, add=add)
  
}



# Integrated tail function using GPD-MLE estimates
# Pi(R)
.IntTailGPD <- function(data, gamma, sigma, R, warnings = TRUE, plot = TRUE, add = FALSE,
                      main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  if (length(R) > 1) {
    stop("R should be a numeric of length 1.")
  }
  
  if ( length(gamma) != length(sigma)) {
    stop("gamma and sigma should have equal length.")
  } 
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  premium[K] <- (K+1) / (n+1) * .IntTailGPD_aux(t=X[n-K], gamma=gamma[K], sigma=sigma[K], endpoint=Inf, R=R)

  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when R >= X[n-K]
  premium[X[n-K]>R] <- NA
  
  if (any(X[n-K]>R) & warnings) {
    warning("R is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R), plot=plot, add=add)
  
}


# Auxiliary function used in ExcessSplice
.IntTailGPD_aux <- function(t, gamma, sigma, R, endpoint = Inf) {
  
  if (any(gamma >= 1)) stop("gamma should be strictly smaller than 1.")
  
  if (is.finite(endpoint)) {
    
    pT <- pgpd(endpoint, mu=t, gamma=gamma, sigma=sigma)
    premium <- ( sigma/(1-gamma) * ( (1 + gamma/sigma * (R-t) )^(1-1/gamma) - (1 + gamma/sigma * (endpoint-t) )^(1-1/gamma) ) + (1-pT)*(R-endpoint)) / pT
    # Special case if R>=endpoint
    premium[R >= endpoint] <- 0
    
  } else {
    premium <- sigma/(1-gamma) * (1 + gamma/sigma * (R-t) )^(1-1/gamma)
  }
  

  return(premium)
}



# Premium of excess-loss insurance with retention R and limit L using Hill estimates (in Pareto model)
# Pi(R) - Pi(R+L)
ExcessPareto <- function(data, gamma, R, L = Inf, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                        main="Estimates for premium of excess-loss insurance", ...) {
  
  if (any(R < 0)) {
    stop("R should be positive.")
  }
  
  if (any(L < 0)) {
    stop("L should be positive.")
  }
  
  if (length(R) != length(L)) {
    if (length(R) != 1 & length(L) != 1) {
      stop("R and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  f <- function(R) .IntTailPareto(R=R, data=data, gamma=gamma, endpoint=endpoint, warnings=warnings, 
                                  plot=FALSE, add=FALSE) 
  
  res <- f(R)
  Im <- res$premium
  K <- res$k
  
  if (is.finite(L)) {
    Il <- f(R+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R, L=L), plot=plot, add=add)
  
}

# Old function name
ExcessHill <- ExcessPareto

# Premium of excess-loss insurance with retention R and limit L using GPD-MLE estimates
# Pi(R) - Pi(R+L)
ExcessGPD <- function(data, gamma, sigma, R, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                       main="Estimates for premium of excess-loss insurance", ...) {
  
  if (any(R < 0)) {
    stop("R should be positive.")
  }
  
  if (any(L < 0)) {
    stop("L should be positive.")
  }
  
  if (length(L) != length(R)) {
    if (length(L) != 1 & length(R) != 1) {
      stop("R and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  f <- function(R) .IntTailGPD(R=R, data=data, gamma=gamma, sigma=sigma, warnings=warnings, 
                               plot=FALSE, add=FALSE) 
  
  res <- f(R)
  Im <- res$premium
  K <- res$k
  
  if (is.finite(L)) {
    Il <- f(R+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R, L=L), plot=plot, add=add)
  
}

# Premium of excess-loss insurance with retention R and limit L using EPD estimates
# Pi(R) - Pi(R+L)
ExcessEPD <- function(data, gamma, kappa, tau, R, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                      main="Estimates for premium of excess-loss insurance", ...) {
  
  
  if (any(R < 0)) {
    stop("R should be positive.")
  }
  
  if (any(L < 0)) {
    stop("L should be positive.")
  }
  
  if (length(L) != length(R)) {
    if (length(L) != 1 & length(R) != 1) {
      stop("R and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  f <- function(R) .IntTailEPD(R=R, data=data, gamma=gamma, kappa=kappa, tau=tau, warnings=warnings, 
                              plot=FALSE, add=FALSE) 
  
  res <- f(R)
  Im <- res$premium
  K <- res$k
  
  if (is.finite(L)) {
    Il <- f(R+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  .plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  .output(list(k=K, premium=premium, R=R, L=L), plot=plot, add=add)
  
}


#######################################################

# Integrated tail function of splicing of ME and (truncated) Pareto
# Pi(R)
.IntTailSplicePareto <- function(R, splicefit) {
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- EVTfit$endpoint
  const <- splicefit$const
  l <- length(const)
  
  trunclower <- splicefit$trunclower
  tvec <- splicefit$t
  
  
  premium <- numeric(length(R))
  
  for (i in 1:length(R)) {
    
    # Logical indicating if final premium is calculated
    end <- FALSE
    
    if (R[i] > tvec[l]) {
      # R[i]>tvec[l] case
      premium[i] <- (1-const[l]) * .IntTailPareto_aux(t=tvec[l], gamma=EVTfit$gamma[l], endpoint=endpoint[l], R=R[i])
      
      end <- TRUE 
      
    } else {
      # Integrate from tvec[l] to endpoint[l]
      premium[i] <- (1-const[l]) * .IntTailPareto_aux(t=tvec[l], gamma=EVTfit$gamma[l], endpoint=endpoint[l], R=tvec[l])
    } 
      
    
    # Only necessary if more than 2 EVT parts in splicing
    if (l>1 & !end) {  
      
      for (j in (l-1):1)  {

        if (R[i]>tvec[j] & R[i] <= tvec[j+1]) {
          # tvec[j]<R[i]<tvec[j+1] case
        
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- .IntTailPareto_aux(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, R=tvec[j+1])
          par_u <- .IntTailPareto_aux(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, R=R[i])
          
          # Integrate from R[i] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_u - par_tt) + (1-const[j+1]) * (tvec[j+1]-R[i]) + premium[i]
          
          
          # Stop for-loop
          break
          
          # Final premium obtained
          end <- TRUE
          
          
        } else {
          
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- .IntTailPareto_aux(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, R=tvec[j+1])
          par_t <- .IntTailPareto_aux(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, R=tvec[j])
          
          # Integrate from tvec[j+1] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_t - par_tt) + (1-const[j+1]) * (tvec[j+1]-tvec[j]) + premium[i]
          
        }
         
      }
      
    }
    
    if (R[i] <= trunclower) {
      premium[i] <- premium[i] + (trunclower-R[i])
      R[i] <- trunclower
    }
    
    if (R[i] <= tvec[1]) {
      # R[i]<tvec[1] case

      # Set C to Inf because C is maximal cover amount
      me_u <- .ME_XL(R=R[i], C=Inf, shape = MEfit$shape, alpha = MEfit$p, 
                     theta = MEfit$theta)
      me_t <- .ME_XL(R=tvec[1], C=Inf, shape = MEfit$shape, alpha = MEfit$p, 
                     theta = MEfit$theta)
      F0 <- .ME_cdf(trunclower, shape = MEfit$shape, alpha = MEfit$p, theta = MEfit$theta)
      Ft <- .ME_cdf(tvec[1], shape = MEfit$shape, alpha = MEfit$p, theta = MEfit$theta)


      # Integrate from R[i] to tvec[1] and add to premium
      # Take truncation at trunclower and tvec[1] into account!
      premium[i] <- ( me_u-me_t + (Ft-1) * (tvec[1]-R[i]) ) / (Ft-F0) * const[1] + (1-const[1]) * (tvec[1]-R[i]) + premium[i]
    }
      
  }
  
  return(premium)
}


# Integrated tail function of splicing of ME and GPD
# Pi(R)
.IntTailSpliceGPD <- function(R, splicefit) {
  
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- EVTfit$endpoint
  const <- splicefit$const
  l <- length(const)
  
  trunclower <- splicefit$trunclower
  tvec <- splicefit$t
  
  premium <- numeric(length(R))
  
  for (i in 1:length(R)) {
    
    # Logical indicating if final premium is calculated
    end <- FALSE
    
    if (R[i] > tvec[l]) {
      # R[i]>tvec[l] case
      premium[i] <- (1-const[l]) * .IntTailGPD_aux(t=tvec[l], gamma=EVTfit$gamma[l], sigma=EVTfit$sigma[l], 
                                                   endpoint=endpoint[l], R=R[i])
      end <- TRUE 
      
    } else {
      # Integrate from tvec[l] to endpoint[l]
      premium[i] <- (1-const[l]) * .IntTailGPD_aux(t=tvec[l], gamma=EVTfit$gamma[l], sigma=EVTfit$sigma[l], 
                                                   endpoint=endpoint[l], R=tvec[l])
    }  
    
    # Only necessary if more than 2 EVT parts in splicing
    if (l>1 & !end) {
      
      for (j in (l-1):1)  {
       
        if (R[i]>tvec[j] & R[i] <= tvec[j+1]) {
		      # tvec[j]<R[i]<tvec[j+1] case 
		
		      # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- .IntTailGPD_aux(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, R=tvec[j+1])
          par_u <- .IntTailGPD_aux(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, R=R[i])
          
          # Integrate from R[i] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_u - par_tt) + (1-const[j+1]) * (tvec[j+1]-R[i]) + premium[i]  
          
          # Stop for-loop
          break
          
          # Final premium obtained
          end <- TRUE
          
          
        } else {
          
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- .IntTailGPD_aux(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, R=tvec[j+1])
          par_t <- .IntTailGPD_aux(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, R=tvec[j])
          
          # Integrate from tvec[j] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_t - par_tt) + (1-const[j+1]) * (tvec[j+1]-tvec[j]) + premium[i]  
        }
        
      }
      
    }
    
    if (R[i] <= trunclower) {
      premium[i] <- premium[i] + (trunclower-R[i])
      R[i] <- trunclower
    }
    
    if (R[i] <= tvec[1]) {
      # R[i]<tvec[1] case
      
      # Set C to Inf because C is maximal cover amount
      me_u <- .ME_XL(R=R[i], C=Inf, shape = MEfit$shape, alpha = MEfit$p, 
                     theta = MEfit$theta)
      me_t <- .ME_XL(R=tvec[1], C=Inf, shape = MEfit$shape, alpha = MEfit$p, 
                     theta = MEfit$theta)
      F0 <- .ME_cdf(trunclower, shape = MEfit$shape, alpha = MEfit$p, theta = MEfit$theta)
      Ft <- .ME_cdf(tvec[1], shape = MEfit$shape, alpha = MEfit$p, theta = MEfit$theta)
      
      # Integrate from R[i] to tvec[1] and add to premium
      # Take truncation at tvec[1] into account!
      premium[i] <- ( me_u-me_t + (Ft-1) * (tvec[1]-R[i]) ) / (Ft-F0) * const[1] + (1-const[1]) * (tvec[1]-R[i]) + premium[i]
    }
    
  }
  
  return(premium)
}


# Premium of excess-loss insurance with retention R and limit L
# Pi(R) - Pi(R+L)
ExcessSplice <- function(R, L=Inf, splicefit) {
  
  if (any(R < 0)) {
    stop("R should be positive.")
  }
  
  if (any(L < 0)) {
    stop("L should be positive.")
  }
  
  if (length(L) != length(R)) {
    if (length(L) != 1 & length(R) != 1) {
      stop("R and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  # Check input
  if (class(splicefit) != "SpliceFit") stop("splicefit should be of class SpliceFit.")
  
  
  type <- splicefit$type
  
  # 2 since type[1]="ME"
  if (type[2] == "GPD") {
    f <- function(R) .IntTailSpliceGPD(R, splicefit=splicefit)
    
  } else if (type[2] %in% c("Pa", "TPa")) {
    f <- function(R) .IntTailSplicePareto(R, splicefit=splicefit)

  } else {
    stop("Invalid type.")
  }
  
  
  
  # First premium
  Im <- f(R)
  
  # Length of premium vector
  le <- max(length(L), length(R))
  
  # Make L vector of length le
  L <- rep(L, length.out=le)
  
  # Second premium
  Il <- f(R+L)
  
  # Numerical issues can arrise when L=Inf
  Il[is.infinite(L)] <- 0
  
  
  # Final premium
  premium <- Im - Il
  
  return(premium)
  
}



#########################################################

# Value-at-risk
VaR <- function(p, splicefit) {
  
  return(qSplice(p=1-p, splicefit=splicefit))
}


# Expected Shortfall (ES)
ES <- function(p, splicefit) {
  
  
  # Check input
  if (!is.numeric(p)) {
    stop("p should be numeric.")
  }
  
  if (any(p <= 0) | any(p > 1)) {
    stop("All elements of p should be in (0,1].")
  }
  
  if (class(splicefit) != "SpliceFit") stop("splicefit should be of class SpliceFit.")
  
  # VaR
  VaR <- VaR(p=p, splicefit=splicefit)

  # ES
  es <-  VaR + 1/p * ExcessSplice(VaR, splicefit=splicefit)
    
  return(es)
}




