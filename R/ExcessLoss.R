

# Premium for excess loss insurance with retention u 
# using Hill estimates
ExcessHill <- function(data, gamma, u, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
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
  if(is.finite(endpoint)) {
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
  
  output(list(k=K, premium=premium, u=u),plot=plot,add=add)
  
}


ExcessHill_single <- function(t, gamma, u, endpoint = Inf) {
  
  # Premium
  if(is.finite(endpoint)) {
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


# Premium for excess loss insurance with retention u 
# using EPD estimates
ExcessEPD <- function(data, gamma, delta, tau, u, warnings = TRUE, plot = TRUE, add = FALSE,
                    main="Estimates for premium of excess-loss insurance", ...) {
  
  # Check input arguments
  checkInput(data)
  
  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  if( length(gamma)!=length(delta) | length(gamma)!=length(tau)) {
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
  
  output(list(k=K, premium=premium, u=u),plot=plot,add=add)
  
}



# Premium for excess loss insurance with retention u 
# using GPD-MLE estimates
ExcessGPD <- function(data, gamma, sigma, u, warnings = TRUE, plot = TRUE, add = FALSE,
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
  
  output(list(k=K, premium=premium, u=u),plot=plot,add=add)
  
}


# Premium for excess loss insurance with retention u 
# using splicing of ME and (truncated) Pareto
ExcessSpliceHill <- function(splicefit, u) {
  

  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit

  endpoint <- Inf
  if(splicefit$type=="trHill") endpoint <- EVTfit$endpoint

  if(exists("const2",where=splicefit)) {
    # Double splicing
    
    if(u<endpoint) {
      
      if(u > splicefit$t2 ) {
        # u>t2 case
        premium <- (1-splicefit$const2) * ExcessHill_single(t=splicefit$t2,gamma=EVTfit$gamma2,endpoint=Inf,u=u)
        
      } else if(u<splicefit$t2 & u>splicefit$t) {
        # t<u<t2 case
        
        par_t2 <- ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=splicefit$t2)
        par_u <- ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=u)
        
        premium <- (splicefit$const2-splicefit$const) * (par_u - par_t2) + (1-splicefit$const2) * (splicefit$t2-u) +  
                    (1-splicefit$const2) * ExcessHill_single(t=splicefit$t2,gamma=EVTfit$gamma2,endpoint=Inf,u=splicefit$t2)
     
      } else {
        # u<t case
        
        # Set C to Inf because C is maximal cover amount
        me_u <- ME_XL(R=u, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                      theta = MEfit$theta)
        me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                      theta = MEfit$theta)
        Ft <- ME_cdf(splicefit$t,shape = MEfit$shape, alpha = MEfit$alpha, theta = MEfit$theta)
        
        par_t <- ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=splicefit$t)
        
        par_t2 <- ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=splicefit$t2)
        
        premium <- (me_u-me_t - (1-Ft) * (splicefit$t-u))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u) + 
          (splicefit$const2-splicefit$const) * (par_t - par_t2) + 
          (1-splicefit$const2) * (splicefit$t2-splicefit$t) + (1-splicefit$const2) * ExcessHill_single(t=splicefit$t2,gamma=EVTfit$gamma2,endpoint=Inf,u=splicefit$t2)
      }
      
    } else {
      premium <- 0
    }
    
    
  } else {
    # Single splicing
    
    if(u > splicefit$t & u<endpoint) {
      # u>t case
      premium <- (1-splicefit$const) * ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=u)
      
    } else if(u<endpoint) {
      # u<t case
      # Set C to Inf because C is maximal cover amount
      me_u <- ME_XL(R=u, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                    theta = MEfit$theta)
      me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                    theta = MEfit$theta)
      Ft <- ME_cdf(splicefit$t,shape = MEfit$shape, alpha = MEfit$alpha, theta = MEfit$theta)
      
      premium <- (me_u-me_t - (1-Ft) * (splicefit$t-u))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u) + 
        (1-splicefit$const) * ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=splicefit$t)
      
    } else {
      premium <- 0
    }
  }
  
 
  return(premium)
}




# Premium for excess loss insurance with retention u 
# using splicing of ME and (truncated) censored Pareto
ExcessSplicecHill <- function(splicefit, u) {
  
  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  endpoint <- Inf
  if(splicefit$type=="trcHill") endpoint <- EVTfit$endpoint
  
  if(u > splicefit$t & u<endpoint) {
    # u>t case
    premium <- (1-splicefit$const) * ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=u)
    
  } else if(u<endpoint) {
    # u<t case
    # Set C to Inf because C is maximal cover amount
    me_u <- ME_XL(R=u, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                  theta = MEfit$theta)
    me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                  theta = MEfit$theta)
    Ft <- ME_cdf(splicefit$t,shape = MEfit$shape, alpha = MEfit$alpha, theta = MEfit$theta)
    
    premium <- (me_u-me_t - (1-Ft) * (splicefit$t-u))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u) + 
      (1-splicefit$const) * ExcessHill_single(t=splicefit$t,gamma=EVTfit$gamma,endpoint=endpoint,u=splicefit$t)
    
  } else {
    premium <- 0
  }
  
  return(premium)
}


# Premium for excess loss insurance with retention u 
# using splicing of ME and GPD
ExcessSpliceGPD <- function(splicefit, u) {
  
  
  if (length(u)>1) {
    stop("u should be a numeric of length 1.")
  }
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit

  if(u > splicefit$t) {
    # u>t case
    premium <- (1-splicefit$const) *  EVTfit$sigma/(1-EVTfit$gamma) * (1 + EVTfit$gamma/EVTfit$sigma * (u-splicefit$t) )^(1-1/EVTfit$gamma)
    
  } else {
    # u<t case
    # Set C to Inf because C is maximal cover amount
    me_u <- ME_XL(R=u, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                  theta = MEfit$theta)
    me_t <- ME_XL(R=splicefit$t, C=Inf, shape = MEfit$shape, alpha = MEfit$alpha, 
                  theta = MEfit$theta)
    Ft <- ME_cdf(splicefit$t,shape = MEfit$shape, alpha = MEfit$alpha, theta = MEfit$theta)

   premium <- (me_u-me_t - (1-Ft) * (splicefit$t-u))/Ft * splicefit$const + (1-splicefit$const) * (splicefit$t-u) + (1-splicefit$const) * EVTfit$sigma / (1-EVTfit$gamma)
      
  } 
  
  return(premium)
}