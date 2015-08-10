



# Integrated tail function using Hill estimates
IntTailHill <- function(data, gamma, u, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                    main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  checkInput(data,gamma)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  
  
  if (length(u)>1 & length(u)!=(n-1)) {
    stop("u should be a numeric of length 1 or n-1.")
  }
  
  # Premium
  if (is.finite(endpoint)) {
    
    if (u>=endpoint) {
      premium <- 0
      
    } else {
      # Truncated Pareto model
      #beta <- X[n]/X[n-K]
      beta <- endpoint/X[n-K]
      premium[K] <- (K+1) / (n+1) / (1-beta^(-1/gamma[K])) * ( X[n-K]^(1/gamma[K]) / (1/gamma[K]-1) * 
                    (u^(1-1/gamma[K]) - endpoint^(1-1/gamma[K])) + beta^(-1/gamma[K]) * (u - endpoint) )
    }
  } else {
    # Pareto model
    premium[K] <- (K+1) / (n+1) / (1/gamma[K]-1) * X[n-K]^(1/gamma[K]) * u^(1-1/gamma[K])
  }

  
  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when u >= X[n-K]
  premium[X[n-K]>u] <- NA
  
  if (any(X[n-K]>u) & warnings) {
    warning("u is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, u=u), plot=plot, add=add)
  
}

# Auxiliary function used in ExcessSplice
IntTailHill_single <- function(t, gamma, u, endpoint = Inf) {
  
  # Premium
  if (is.finite(endpoint)) {
    
    if (u>=endpoint) {
      premium <- 0
      
    } else {
      # Truncated Pareto model
      beta <- endpoint/t
      premium <- 1 / (1-beta^(-1/gamma)) * ( t^(1/gamma) / (1/gamma-1) * 
                                               (u^(1-1/gamma) - endpoint^(1-1/gamma)) + beta^(-1/gamma) * (u - endpoint) )
    }
    
  } else {
    # Pareto model
    premium <- 1 / (1/gamma-1) * t^(1/gamma) * u^(1-1/gamma)
  }
  return(premium)
}


# Integrated tail function using EPD estimates
IntTailEPD <- function(data, gamma, delta, tau, u, warnings = TRUE, plot = TRUE, add = FALSE,
                    main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  checkInput(data)
  
  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  if (length(gamma)!=length(delta) | length(gamma)!=length(tau)) {
    stop("gamma, delta and tau should have equal length.")
  } 
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  premium[K] <- (K+1) / (n+1) * X[n-K]^(1/gamma[K]) * ( (1-delta[K]/gamma[K]) / (1/gamma[K]-1) * u^(1-1/gamma[K]) +
                                                        delta[K]/(gamma[K]*X[n-K]^tau[K]) / (1/gamma[K]-tau[K]-1) *
                                                          u^(1+tau[K]-1/gamma[K]) )
  
  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when u >= X[n-K]
  premium[X[n-K]>u] <- NA
  
  if (any(X[n-K]>u) & warnings) {
    warning("u is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, u=u), plot=plot, add=add)
  
}



# Integrated tail function using GPD-MLE estimates
IntTailGPD <- function(data, gamma, sigma, u, warnings = TRUE, plot = TRUE, add = FALSE,
                      main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  checkInput(data)
  
  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  if( length(gamma)!=length(sigma)) {
    stop("gamma and sigma should have equal length.")
  } 
  
  
  X <- as.numeric(sort(data))
  n <- length(X)
  premium <- numeric(n)
  K <- 1:(n-1)
  
  premium[K] <- (K+1) / (n+1) * sigma[K]/(1-gamma[K]) * (1 + gamma[K]/sigma[K] * (u-X[n-K]) )^(1-1/gamma[K])

  # Remove negative values
  premium[premium<0] <- NA
  
  # Only valid when u >= X[n-K]
  premium[X[n-K]>u] <- NA
  
  if (any(X[n-K]>u) & warnings) {
    warning("u is smaller than X[n-k] for some K, use global fits for these cases!")
  }
  
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, u=u), plot=plot, add=add)
  
}


# Auxiliary function used in ExcessSplice
IntTailGPD_single <- function(t, gamma, sigma, u, endpoint = Inf) {
  
  if (is.finite(endpoint)) {
    
    if (u>=endpoint) {
      premium <- 0
      
    } else {
      pT <- pgpd(endpoint, mu=t, gamma=gamma, sigma=sigma)
      premium <- ( sigma/(1-gamma) * ( (1 + gamma/sigma * (u-t) )^(1-1/gamma) - (1 + gamma/sigma * (endpoint-t) )^(1-1/gamma) ) + (1-pT)*(u-endpoint)) / pT
    }
    
  } else {
    premium <- sigma/(1-gamma) * (1 + gamma/sigma * (u-t) )^(1-1/gamma)
  }
  

  return(premium)
}



# Premium of excess-loss insurance with retention M and limit L using Hill estimates
ExcessHill <- function(data, gamma, M, L = Inf, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                        main="Estimates for premium of excess-loss insurance", ...) {
  
  if (any(M<0)) {
    stop("M should be positive.")
  }
  
  if (any(L<0)) {
    stop("L should be positive.")
  }
  
  if (length(M)!=length(L)) {
    if(length(M)!=1 & length(L)!=1) {
      stop("M and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  f <- function(u) IntTailHill(u=u, data=data, gamma=gamma, endpoint=endpoint, warnings=warnings, 
                               plot=FALSE, add=FALSE) 
  
  res <- f(M)
  Im <- res$premium
  K <- res$k
  
  if(is.finite(L)) {
    Il <- f(M+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, M=M, L=L), plot=plot, add=add)
  
}


# Premium of excess-loss insurance with retention M and limit L using GPD-MLE estimates
ExcessGPD <- function(data, gamma, sigma, M, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                       main="Estimates for premium of excess-loss insurance", ...) {
  
  if (any(M<0)) {
    stop("M should be positive.")
  }
  
  if (any(L<0)) {
    stop("L should be positive.")
  }
  
  f <- function(u) IntTailGPD(u=u, data=data, gamma=gamma, sigma=sigma, warnings=warnings, 
                               plot=FALSE, add=FALSE) 
  
  res <- f(M)
  Im <- res$premium
  K <- res$k
  
  if(is.finite(L)) {
    Il <- f(M+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, M=M, L=L), plot=plot, add=add)
  
}

# Premium of excess-loss insurance with retention M and limit L using EPD estimates
ExcessEPD <- function(data, gamma, delta, tau, M, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                      main="Estimates for premium of excess-loss insurance", ...) {
  
  
  if (any(M<0)) {
    stop("M should be positive.")
  }
  
  if (any(L<0)) {
    stop("L should be positive.")
  }
  
  f <- function(u) IntTailEPD(u=u, data=data, gamma=gamma, delta=delta, tau=tau, warnings=warnings, 
                              plot=FALSE, add=FALSE) 
  
  res <- f(M)
  Im <- res$premium
  K <- res$k
  
  if(is.finite(L)) {
    Il <- f(M+L)$premium
  } else {
    Il <- 0
  }
  
  premium <- Im - Il
  
  
  # plots if TRUE 
  plotfun(K, premium[K], type="l", xlab="k", ylab="Premium", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, premium=premium, M=M, L=L), plot=plot, add=add)
  
}


#######################################################

# Integrated tail function of splicing of ME and (truncated) Pareto
IntTailSpliceHill <- function(u, splicefit) {
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- EVTfit$endpoint
  const <- splicefit$const
  l <- length(const)
  
  tvec <- splicefit$t
  
  
  premium <- numeric(length(u))
  
  for (i in 1:length(u)) {
    
    # Logical indicating if final premium is calculated
    end <- FALSE
    
    if (u[i] > tvec[l]) {
      # u[i]>tvec[l] case
      premium[i] <- (1-const[l]) * IntTailHill_single(t=tvec[l], gamma=EVTfit$gamma[l], endpoint=endpoint[l], u=u[i])
      
      end <- TRUE 
      
    } else {
      # Integrate from tvec[l] to endpoint[l]
      premium[i] <- (1-const[l]) * IntTailHill_single(t=tvec[l], gamma=EVTfit$gamma[l], endpoint=endpoint[l], u=tvec[l])
    } 
      
    
    # Only necessary if more than 2 EVT parts in splicing
    if (l>1 & !end) {  
      
      for (j in (l-1):1)  {

        if (u[i]>tvec[j] & u[i]<=tvec[j+1]) {
          # tvec[j]<u[i]<tvec[j+1] case
        
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- IntTailHill_single(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, u=tvec[j+1])
          par_u <- IntTailHill_single(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, u=u[i])
          
          # Integrate from u[i] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_u - par_tt) + (1-const[j+1]) * (tvec[j+1]-u[i]) + premium[i]
          
          
          # Stop for-loop
          break
          
          # Final premium obtained
          end <- TRUE
          
          
        } else {
          
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- IntTailHill_single(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, u=tvec[j+1])
          par_t <- IntTailHill_single(t=tvec[j], gamma=EVTfit$gamma[j], endpoint=e, u=tvec[j])
          
          # Integrate from tvec[j+1] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_t - par_tt) + (1-const[j+1]) * (tvec[j+1]-tvec[j]) + premium[i]
          
        }
         
      }
      
    }
          
    if (u[i] <= tvec[1]) {
      # u[i]<tvec[1] case

      # Set C to Inf because C is maximal cover amount
      me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      me_t <- ME_XL(R=tvec[1], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      Ft <- ME_cdf(tvec[1], shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
      
      # Integrate from u[i] to tvec[1] and add to premium
      # Take truncation at tvec[1] into account!
      premium[i] <- (me_u-me_t - (1-Ft) * (tvec[1]-u[i]))/Ft * const[1] + (1-const[1]) * (tvec[1]-u[i]) + premium[i]
    }
      
  }
  
  return(premium)
}


# Integrated tail function of splicing of ME and GPD
IntTailSpliceGPD <- function(u, splicefit) {
  
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- EVTfit$endpoint
  const <- splicefit$const
  l <- length(const)
  
  tvec <- splicefit$t
  
  premium <- numeric(length(u))
  
  for (i in 1:length(u)) {
    
    # Logical indicating if final premium is calculated
    end <- FALSE
    
    if (u[i] > tvec[l]) {
      # u[i]>tvec[l] case
      premium[i] <- (1-const[l]) * IntTailGPD_single(t=tvec[l], gamma=EVTfit$gamma[l], sigma=EVTfit$sigma[l], 
                                                     endpoint=endpoint[l], u=u[i])
      end <- TRUE 
      
    } else {
      # Integrate from tvec[l] to endpoint[l]
      premium[i] <-  (1-const[l]) * IntTailGPD_single(t=tvec[l], gamma=EVTfit$gamma[l], sigma=EVTfit$sigma[l], 
                                                      endpoint=endpoint[l], u=tvec[l])
    }  
    
    # Only necessary if more than 2 EVT parts in splicing
    if (l>1 & !end) {
      
      for (j in (l-1):1)  {
       
        if (u[i]>tvec[j] & u[i]<=tvec[j+1]) {
		      # tvec[j]<u[i]<tvec[j+1] case 
		
		      # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- IntTailGPD_single(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, u=tvec[j+1])
          par_u <- IntTailGPD_single(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, u=u[i])
          
          # Integrate from u[i] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_u - par_tt) + (1-const[j+1]) * (tvec[j+1]-u[i]) + premium[i]  
          
          # Stop for-loop
          break
          
          # Final premium obtained
          end <- TRUE
          
          
        } else {
          
          # Actual endpoint of splicing part
          e <- min(tvec[j+1], endpoint[j])
          
          par_tt <- IntTailGPD_single(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, u=tvec[j+1])
          par_t <- IntTailGPD_single(t=tvec[j], gamma=EVTfit$gamma[j], sigma=EVTfit$sigma[j], endpoint=e, u=tvec[j])
          
          # Integrate from tvec[j] to tvec[j+1] and add to premium
          premium[i] <- (const[j+1]-const[j]) * (par_t - par_tt) + (1-const[j+1]) * (tvec[j+1]-tvec[j]) + premium[i]  
        }
        
      }
      
    }
    
    if (u[i] <= tvec[1]) {
      # u[i]<tvec[1] case
      
      # Set C to Inf because C is maximal cover amount
      me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      me_t <- ME_XL(R=tvec[1], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      Ft <- ME_cdf(tvec[1], shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
      
      # Integrate from u[i] to tvec[1] and add to premium
      # Take truncation at tvec[1] into account!
      premium[i] <- (me_u-me_t - (1-Ft) * (tvec[1]-u[i]))/Ft * const[1] + (1-const[1]) * (tvec[1]-u[i]) + premium[i]
    }
    
  }
  
  return(premium)
}


# Premium of excess-loss insurance with retention M and limit L
ExcessSplice <- function(M, L=Inf, splicefit) {
  
  if (any(M<0)) {
    stop("M should be positive.")
  }
  
  if (any(L<0)) {
    stop("L should be positive.")
  }
  
  if (length(L)!=length(M)) {
    if(length(L)!=1 & length(M)!=1) {
      stop("M and L should have equal length or at least one of them should have length 1.")
    }
  }
  
  
  type <- splicefit$type
  
  # 2 since type[1]="ME"
  if (type[2]=="GPD") {
    f <- function(u) IntTailSpliceGPD(u, splicefit=splicefit)
    
  } else if (type[2] %in% c("Pa", "tPa", "cPa", "ciPa", "trciPa")) {
    f <- function(u) IntTailSpliceHill(u, splicefit=splicefit)

  } else {
    stop("Invalid type.")
  }
  
  
  Im <- f(M)
  
  Il <- f(M+L)
  Il[!is.finite(L)] <- 0
  
  premium <- Im - Il
  
  return(premium)
  
}



#########################################################

# Value-at-risk
VaR <- function(p, splicefit) {
  
  return(SpliceQuant(p=1-p, splicefit=splicefit))
}


# Expected Shortfall (ES)
ES <- function(p, splicefit) {
  
  # Check input
  if (!is.numeric(p)) {
    stop("p should be numeric.")
  }
  
  if (any(p<=0) | any(p>1)) {
    stop("All elements of p should be in (0,1].")
  }
  
  # VaR
  VaR <- VaR(p=p, splicefit=splicefit)

  # Conditional tail expectation 
  es <-  VaR + 1/p * ExcessSplice(VaR, splicefit=splicefit)
    
  return(es)
}




