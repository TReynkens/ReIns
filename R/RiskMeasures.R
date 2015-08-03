



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
    # Truncated Pareto model
    #beta <- X[n]/X[n-K]
    beta <- endpoint/X[n-K]
    premium[K] <- (K+1) / (n+1) / (1-beta^(-1/gamma[K])) * ( X[n-K]^(1/gamma[K]) / (1/gamma[K]-1) * 
                  (u^(1-1/gamma[K]) - endpoint^(1-1/gamma[K])) + beta^(-1/gamma[K]) * (u - endpoint) )
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
    # Truncated Pareto model
    beta <- endpoint/t
    premium <- 1 / (1-beta^(-1/gamma)) * ( t^(1/gamma) / (1/gamma-1) * 
                                                               (u^(1-1/gamma) - endpoint^(1-1/gamma)) + beta^(-1/gamma) * (u - endpoint) )
  } else {
    # Pareto model
    premium <-1 / (1/gamma-1) * t^(1/gamma) * u^(1-1/gamma)
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


# Premium of excess-loss insurance with retentation M and limit L using Hill estimates
ExcessHill <- function(data, gamma, M, L = Inf, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                        main="Estimates for premium of excess-loss insurance", ...) {
  
  if (M<0) {
    stop("M should be positive.")
  }
  
  if (L<0) {
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


# Premium of excess-loss insurance with retentation M and limit L using GPD-MLE estimates
ExcessGPD <- function(data, gamma, sigma, M, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                       main="Estimates for premium of excess-loss insurance", ...) {
  
  if (M<0) {
    stop("M should be positive.")
  }
  
  if (L<0) {
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

# Premium of excess-loss insurance with retentation M and limit L using EPD estimates
ExcessEPD <- function(data, gamma, delta, tau, M, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
                      main="Estimates for premium of excess-loss insurance", ...) {
  
  
  if (M<0) {
    stop("M should be positive.")
  }
  
  if (L<0) {
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
  
  if(l>2) {
    stop("Splicing of ME with 3 or more Pareto distributions is not supported.")
  }
  
  premium <- numeric(length(u))
  
  for (i in 1:length(u)) {
    
    if (l>1) {
      # l-splicing
      
      if (u[i]<endpoint[2]) {
        
        if (u[i] > tvec[2]) {
          # u[i]>t2 case
          premium[i] <- (1-const[2]) * IntTailHill_single(t=tvec[2], gamma=EVTfit$gamma[2], endpoint=endpoint[2], u=u[i])
          
        } else if (u[i]<tvec[2] & u[i]>tvec[1]) {
          # t<u[i]<t2 case
          
          e <- min(tvec[2],endpoint[1])
          
          par_t2 <- IntTailHill_single(t=tvec[1], gamma=EVTfit$gamma[1], endpoint=e, u=tvec[2])
          par_u <- IntTailHill_single(t=tvec[1], gamma=EVTfit$gamma[1], endpoint=e, u=u[i])
          
          premium[i] <- (const[2]-const[1]) * (par_u - par_t2) + (1-const[2]) * (tvec[2]-u[i]) +  
            (1-const[2]) * IntTailHill_single(t=tvec[2], gamma=EVTfit$gamma[2], endpoint=endpoint[2], u=tvec[2])
          
        } else {
          # u<t case
          
          # Set C to Inf because C is maximal cover amount
          me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                        theta = MEfit$theta)
          me_t <- ME_XL(R=tvec[1], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                        theta = MEfit$theta)
          Ft <- ME_cdf(tvec[1], shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
          
          e <- min(tvec[2],endpoint[1])
          par_t <- IntTailHill_single(t=tvec[1], gamma=EVTfit$gamma[1], endpoint=e, u=tvec[1])
          par_t2 <- IntTailHill_single(t=tvec[1], gamma=EVTfit$gamma[1], endpoint=e, u=tvec[2])
          
          premium[i] <- (me_u-me_t - (1-Ft) * (tvec[1]-u[i]))/Ft * const[1] + (1-const[1]) * (tvec[1]-u[i]) + 
            (const[2]-const[1]) * (par_t - par_t2) + 
            (1-const[2]) * (tvec[2]-tvec[1]) + 
            (1-const[2]) * IntTailHill_single(t=tvec[2], gamma=EVTfit$gamma[2], endpoint=endpoint[2], u=tvec[2])
        }
        
      } else {
        premium[i] <- 0
      }
      
      
    } else {
      # Single splicing
      
      if (u[i]>tvec[l] & u[i]<endpoint[l]) {
        # u[i]>t case
        premium[i] <- (1-const[l]) * IntTailHill_single(t=tvec[l], gamma=EVTfit$gamma, endpoint=endpoint, u=u[i])
        
      } else if (u[i]<endpoint[l]) {
        # u[i]<t case
        # Set C to Inf because C is maximal cover amount
        me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                      theta = MEfit$theta)
        me_t <- ME_XL(R=tvec[l], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                      theta = MEfit$theta)
        Ft <- ME_cdf(tvec[l], shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
        
        premium[i] <- (me_u-me_t - (1-Ft) * (tvec[l]-u[i]))/Ft * const[l] + (1-const[l]) * (tvec[l]-u[i]) + 
          (1-const[l]) * IntTailHill_single(t=tvec[l], gamma=EVTfit$gamma, endpoint=endpoint[l], u=tvec[l])
        
      } else {
        premium[i] <- 0
      }
    }
    
  }
  
  
  return(premium)
}



# Integrated tail function of splicing of ME and (truncated) censored Pareto
IntTailSplicecHill <- function(u, splicefit) {
  

  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- Inf
  if(splicefit$type=="trcHill") endpoint <- EVTfit$endpoint
  
  premium <- numeric(length(u))
  
  for (i in 1:length(u)) {
    
    if (u[i]>splicefit$t & u[i]<endpoint) {
      # u>t case
      premium[i] <- (1-splicefit$const) * IntTailHill_single(t=splicefit$t, gamma=EVTfit$gamma, endpoint=endpoint, u=u[i])
      
    } else if (u[i]<endpoint) {
      # u<t case
      # Set C to Inf because C is maximal cover amount
      me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      Ft <- ME_cdf(splicefit$t, shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
      
      premium[i] <- (me_u-me_t - (1-Ft) * (splicefit$t-u[i]))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u[i]) + 
        (1-splicefit$const) * IntTailHill_single(t=splicefit$t, gamma=EVTfit$gamma, endpoint=endpoint, u=splicefit$t)
      
    } else {
      premium[i] <- 0
    }
  }
  return(premium)
}


# Integrated tail function of splicing of ME and GPD
IntTailSpliceGPD <- function(u, splicefit) {
  
  

  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit

  premium <- numeric(length(u))
  
  for (i in 1:length(u)) {
    
    if(u[i] > splicefit$t) {
      # u>t case
      premium[i] <- (1-splicefit$const) *  EVTfit$sigma/(1-EVTfit$gamma) * (1 + EVTfit$gamma/EVTfit$sigma * (u[i]-splicefit$t) )^(1-1/EVTfit$gamma)
      
    } else {
      # u<t case
      # Set C to Inf because C is maximal cover amount
      me_u <- ME_XL(R=u[i], C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$beta, 
                    theta = MEfit$theta)
      Ft <- ME_cdf(splicefit$t, shape = MEfit$shape, alpha = MEfit$beta, theta = MEfit$theta)
  
     premium[i] <- (me_u-me_t - (1-Ft) * (splicefit$t-u[i]))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u[i]) + 
       (1-splicefit$const) * EVTfit$sigma / (1-EVTfit$gamma)
        
    } 
  
  }
  return(premium)
}


# Premium of excess-loss insurance with retention M and limit L
ExcessSplice <- function(M, L=Inf, splicefit) {
  
  if (M<0) {
    stop("M should be positive.")
  }
  
  if (L<0) {
    stop("L should be positive.")
  }
  
  
  type <- splicefit$type
  
  if(type[1]=="GPD") {
    f <- function(u) IntTailSpliceGPD(u, splicefit=splicefit)
    
  } else if(type[1] %in% c("Hill", "trHill")) {
    f <- function(u) IntTailSpliceHill(u, splicefit=splicefit)
    
  } else if(type[1] %in% c("cHill", "ciHill", "trciHill")) {
    f <- function(u) IntTailSplicecHill(u, splicefit=splicefit)
    
  } else {
    stop("Invalid type.")
  }
  
  
  Im <- f(M)
  
  if(is.finite(L)) {
    Il <- f(M+L)
  } else {
    Il <- 0
  }
  
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




