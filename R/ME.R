######################################################################
## EM algorithm mixtures of Erlangs for censored and truncated data ##
######################################################################

# Author: Roel Verbelen with adaptations by Tom Reynkens

## Initial values

# supply number of shapes M and spread factor s

.ME_initial <- function(lower, upper, trunclower = 0, truncupper = Inf, M = 10, s = 1) {
  
  # data for initial step: treat left and right censored data as observed and take mean for interval censored data
  uc <- lower[lower==upper & !is.na(lower==upper)]
  lc <- upper[is.na(lower)]
  rc <- lower[is.na(upper)]
  ic <- (lower[lower!=upper & !is.na(lower!=upper)] + upper[lower!=upper & !is.na(lower!=upper)]) / 2
  initial_data = c(uc, lc, rc, ic)
  
  # Initial value of theta using spread factor s
  theta <- max(initial_data, na.rm=TRUE) / s
  
  # Initial shapes as quantiles
  shape <- unique(ceiling(quantile(initial_data, probs = seq(0, 1, length.out = M), na.rm = TRUE)/theta))
  
  # Initial value for alpha's
  alpha <- rep(0,length(shape))
  alpha[1] <- sum(initial_data <= shape[1]*theta)
  if (length(shape)>1) {
    for (i in 2:length(shape)){
      alpha[i] <- sum(initial_data <= shape[i]*theta & initial_data > shape[i-1]*theta)
    }
  }
  
  # Keep strictly positive alpha's and corresponding shapes
  shape <- shape[alpha>0]
  alpha <- alpha[alpha>0]/sum(alpha)
  
  # alpha to beta
  t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)  
  beta = alpha * t_probabilities / sum(alpha*t_probabilities)
  list(theta=theta, shape=shape, alpha=alpha, beta=beta) 
}

## Log likelihood

.ME_loglikelihood <- function(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring) {
  
  likelihood_contribution <- numeric(0)
  if(no_censoring){  
    # matrix containing alpha*density (uncensored)
    x_components <-  t(t(x_densities)*beta/t_probabilities)
    # likelihood contribution (uncensored)
    likelihood_contribution <- rowSums(x_components) 
  }   
  if(censoring){  
    # matrix containing alpha*probabilities (censored)
    c_components <- t(t(c_probabilities)*beta/t_probabilities)
    # likelihood contribution (censored)
    likelihood_contribution <- c(likelihood_contribution, rowSums(c_components)) 
  }   
  loglikelihood_contribution <- ifelse(likelihood_contribution>0, log(likelihood_contribution), -1000)
  # loglikelihood
  sum(loglikelihood_contribution)
}

## ^{u}z_{ij}^{(k)}: posterior probabilities (uncensored)

.ME_u_z <-function(x_densities, beta, t_probabilities, M) {      
  
  x_components <- t(t(x_densities)*beta/t_probabilities)
  u_z <- x_components / rowSums(x_components)
  # in case all ^{u}z_{ij}^{k} for j=1,...,M are numerically 0
  u_z[is.nan(u_z)] = 1/M
  u_z
}

## ^{c}z_{ij}^{(k)}: posterior probabilities (censored)

.ME_c_z <-function(c_probabilities, beta, t_probabilities, M) { 
  
  c_components <- t(t(c_probabilities)*beta/t_probabilities)
  c_z <- c_components / rowSums(c_components)
  # in case all ^{c}z_{ij}^{k} for j=1,...,M are numerically 0
  c_z[is.nan(c_z)] = 1/M
  c_z
}
 
## Expected value of censored observations

.ME_expected_c <-function(lower, upper, shape, theta, c_z) {  
  
  c_exp <- theta * (outer(upper, shape+1, pgamma, scale=theta) - outer(lower, shape+1, pgamma, scale=theta))/(outer(upper, shape, pgamma, scale=theta) - outer(lower, shape, pgamma, scale=theta))
  c_exp <- t(t(c_exp)*shape)
  # replace numerical 0/0 (NaN) or Inf by correct expected value
  c_exp <- ifelse(is.nan(c_exp) | c_exp==Inf, ifelse(outer(lower, shape*theta, ">"), lower, upper), c_exp)
  c_exp <- rowSums(c_z * c_exp)  
}

## T^{(k)}: correction term due to truncation

.ME_T <- function(trunclower, truncupper, shape, theta, beta) {
  
  # avoid NaN
  if(truncupper==Inf){ 
    # take log first for numerical stability (avoid Inf / Inf)
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(1 - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)    
  } else{    
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc_log_2 <- shape*log(truncupper)-truncupper/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)-exp(deriv_trunc_log_2)
  }
  sum(beta*deriv_trunc)
}

## Auxiliary functions used to maximize theta in the M-step

# This function returns log of optimal theta!
.theta_nlm <- function(ltheta, x, c_exp, n, beta, shape, trunclower, truncupper) {
  # Ensure theta is strictly positive
  theta <- exp(ltheta)
  TT <- .ME_T(trunclower, truncupper, shape, theta, beta)
  (theta - ((sum(x)+sum(c_exp))/n-TT)/sum(beta*shape))^2
}



## EM algorithm (censored and truncated data)

.ME_em <- function(lower, upper, trunclower = 0, truncupper = Inf, theta, shape, beta, eps = 10^(-3), 
                   beta_tol = 10^(-5), maxiter = Inf) {
  
  n <- length(lower)
  M <- length(shape)
  # separate uncensored and censored observations
  uncensored <- (lower==upper & !is.na(lower==upper))
  # Uncensored observations
  x <- lower[uncensored]
  # Censored observations
  lower[is.na(lower)] <- trunclower
  upper[is.na(upper)] <- truncupper
  lower <- lower[!uncensored]
  upper <- upper[!uncensored]  
  # Boolean for having uncensored and censored observations
  no_censoring <- (length(x) != 0)  
  censoring <- (length(lower) != 0)  
  iteration <- 1
  if(no_censoring){  
    # matrix containing densities (uncensored)
    x_densities <- outer(x,shape,dgamma, scale=theta)
  }   
  if(censoring){  
    # matrix containing censoring probabilities (censored)
    c_probabilities <- outer(upper,shape,pgamma, scale=theta)-outer(lower,shape,pgamma, scale=theta)
  } 
  # truncation probabilities
  t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
  loglikelihood <- .ME_loglikelihood(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring)
  old_loglikelihood <- -Inf
  
  while(loglikelihood - old_loglikelihood > eps & iteration <= maxiter){
    old_loglikelihood <- loglikelihood
    # E step
    if(no_censoring & censoring){
      u_z <- .ME_u_z(x_densities, beta, t_probabilities, M)
      c_z <- .ME_c_z(c_probabilities, beta, t_probabilities, M)
      c_exp <- .ME_expected_c(lower, upper, shape, theta, c_z)      
    } else if(no_censoring){
      u_z <- .ME_u_z(x_densities, beta, t_probabilities, M)
      c_z <- as.matrix(c(0,0))
      c_exp <- as.matrix(c(0,0))
    } else{
      u_z <- as.matrix(c(0,0))
      c_z <- .ME_c_z(c_probabilities, beta, t_probabilities, M)
      c_exp <- .ME_expected_c(lower, upper, shape, theta, c_z)
    }
    # M step
    
    beta <- (colSums(u_z)+colSums(c_z))/n 

    # Remove small beta's
    if(min(beta) < beta_tol){
      shape <- shape[beta > beta_tol]
      beta <- beta[beta > beta_tol]
      beta <- beta/sum(beta)
    }
    
    theta <- exp(nlm(.theta_nlm, log(theta), x, c_exp, n, beta, shape, trunclower, truncupper)$estimate)

    iteration <- iteration + 1
    if(no_censoring){  
      # matrix containing densities (uncensored)
      x_densities <- outer(x,shape,dgamma, scale=theta)
    }   
    if(censoring){  
      # matrix containing censoring probabilities (censored)
      c_probabilities <- outer(upper,shape,pgamma, scale=theta)-outer(lower,shape,pgamma, scale=theta)
    } 
     # truncation probabilities
    t_probabilities <- pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
    loglikelihood <- .ME_loglikelihood(x_densities, c_probabilities, beta, t_probabilities, no_censoring, censoring)
    
  }
  # beta to alpha
  alpha_tilde <- beta / t_probabilities
  alpha <- alpha_tilde / sum(alpha_tilde)
 list(alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, iteration = iteration, AIC=-2*loglikelihood+2*(2*length(alpha)+1),BIC=-2*loglikelihood+(2*length(alpha)+1)*log(n)) 
}

## Shape adjustments

.ME_shape_adj <- function(lower, upper, trunclower = 0, truncupper = Inf, theta, shape, beta, 
                          eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
  shape <- shape[beta>0]
  beta <- beta[beta>0]/sum(beta[beta>0])
  M <- length(shape)
  fit <- .ME_em(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                theta=theta, shape=shape, beta=beta,
                eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  loglikelihood <- fit$loglikelihood
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta
  alpha <- fit$alpha
  M <- length(beta)
  
  # before and after are the loglikelihoods used in the outer while loop
  before_loglikelihood <- -Inf
  after_loglikelihood <- loglikelihood    
  iteration <- 1
  while(after_loglikelihood > before_loglikelihood + eps){    
    
    before_loglikelihood <- after_loglikelihood
    # Try increasing the shapes
    for(i in M:1){
      improve <- TRUE
      while( (improve==TRUE) && (i == M || ifelse(i<=length(shape), shape[i] < shape[i+1]-1, FALSE))) {
        new_shape <- shape
        new_shape[i] <- new_shape[i]+1        
        fit <- .ME_em(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                      theta=theta, shape=new_shape, beta=beta,
                      eps=eps, beta_tol=beta_tol, maxiter=maxiter)
        new_loglikelihood <- fit$loglikelihood
        if(new_loglikelihood > loglikelihood + eps){
          loglikelihood <- new_loglikelihood
          shape <- fit$shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          # number of shape combinations might have changed after EM algorithm if beta_tol > 0
          M <- length(shape)
          i <- min(i, M)
          
        } else {
          improve <- FALSE
        }
      }
    }
    # Try decreasing the shapes
    for(i in 1:M){
      improve <- TRUE
      while( (improve==TRUE) && ( (i == 1) || ifelse(i<=length(shape), shape[i] > shape[i-1]+1, FALSE) ) && ifelse(i<=length(shape), shape[i]>1, TRUE)){
        new_shape <- shape
        new_shape[i] <- new_shape[i]-1
        fit <- .ME_em(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                      theta=theta, shape=new_shape, beta=beta,
                      eps=eps, beta_tol=beta_tol, maxiter=maxiter)
        new_loglikelihood <- fit$loglikelihood
        if(new_loglikelihood > loglikelihood + eps){
          loglikelihood <- new_loglikelihood
          shape <- fit$shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          # number of shape combinations might have changed after EM algorithm if beta_tol > 0
          M <- length(shape)
          i <- min(i, M)
          
        } else { 
          improve <- FALSE
        }
      }          
    }
    after_loglikelihood <- loglikelihood  
    iteration <- iteration + 1
  }
  list(alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(2*length(alpha)+1),BIC=-2*loglikelihood+(2*length(alpha)+1)*log(length(lower))) 
}

## Reduction of M based on an information criterium: AIC and BIC implemented

.ME_shape_red <- function(lower, upper, trunclower = 0, truncupper = Inf, theta, shape, beta, criterium = "AIC", 
                          improve = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, adj = TRUE) {
  n <- length(lower)
  if (adj) {
    fit <- .ME_shape_adj(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                         theta=theta, shape=shape, beta=beta, eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  } else {
    fit <- .ME_em(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                  theta=theta, shape=shape, beta=beta, eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  }
  loglikelihood <- fit$loglikelihood
  IC <- fit[[criterium]]
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta   
  alpha <- fit$alpha
  M <- length(shape)

  while((improve==TRUE) && length(shape) > 1){    
    new_shape <- shape[beta != min(beta)]
    new_beta <- beta[beta != min(beta)]
    new_beta <- new_beta/sum(new_beta)
    if (adj) {
      fit <- .ME_shape_adj(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                           theta=theta, shape=new_shape, beta=new_beta,
                           eps=eps, beta_tol=beta_tol, maxiter=maxiter)
    } else {
      fit <- .ME_em(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                    theta=theta, shape=new_shape, beta=new_beta,
                    eps=eps, beta_tol=beta_tol, maxiter=maxiter)
    }
    new_IC <- fit[[criterium]]
    if(new_IC < IC){ 
      IC <- new_IC
      loglikelihood <- fit$loglikelihood  
      shape <- fit$shape
      theta <- fit$theta
      beta <- fit$beta  
      alpha <- fit$alpha
      M <- length(shape)
     
    } else { 
      improve <- FALSE
    }        
  }
  list(M = M, alpha = alpha, beta = beta, shape = shape, theta = theta, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(2*length(alpha)+1),BIC=-2*loglikelihood+(2*length(alpha)+1)*log(n)) 
}

## Calibration procedure for mixtures of Erlangs by repeatedly using the EM algorithm while adjusting and reducing the shape parameters based on an information criterium (AIC and BIC implemented)

# Specify lower and upper censoring points (lower, upper), lower and upper truncation points (trunclower, truncupper).
# The censoring status is determined as follows:
# Uncensored: lower and upper are both present and equal.
# Left Censored: lower is missing (NA), but upper is present.
# Right Censored: lower is present, but upper is missing (NA).
# Interval Censored: lower and upper are present and different.
# e.g.: lower=c(1,NA,3,4); upper=c(1,2,NA,5); specifies an observed event at 1, left censoring at 2, right censoring at 3, and interval censoring at [4,5],
# By default no truncation: trunclower=0, truncupper=Inf
# alpha = beta in case of no truncation

.ME_fit <- function(lower, upper = lower, trunclower = 0, truncupper = Inf, M = 10, s = 1,
                    criterium = "AIC", reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
  
  initial <- .ME_initial(lower, upper, trunclower, truncupper, M, s)
  # Reduction of the shape combinations
  fit <- .ME_shape_red(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                       theta=initial$theta, shape=initial$shape, beta=initial$beta, 
                       criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=FALSE)
  # Subsequent adjustment and reduction of the shape combinations
  fit <- .ME_shape_red(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                       theta=fit$theta, shape=fit$shape, beta=fit$beta, 
                       criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=TRUE)
  list(alpha = fit$alpha, beta = fit$beta, shape = fit$shape, theta = fit$theta, loglikelihood = fit$loglikelihood, AIC=fit$AIC, BIC=fit$BIC, M = fit$M, M_initial = M, s = s) 
}



## Check input arguments for .MEtune
.ME_checkInput  <- function(lower, upper, trunclower, truncupper, eps, beta_tol, maxiter) {

  nl <- length(lower)
  nu <- length(upper)
  ntl <- length(trunclower)
  ntu <- length(truncupper)
  
  # Check lengths
  if(nl!=1 & nu!=1 & nl!=nu) {
    stop("lower and upper should have equal length if both do not have length 1.")
  }
  
  if(ntl!=1 & ntu!=1 & ntl!=ntu) {
    stop("trunclower and truncupper should have equal length if both do not have length 1.")
  }
  
  
  
  # Check data types
  if(any(!is.numeric(lower) & !is.na(lower))) {
    stop("lower should consist of numerics and/or NAs.")
  }
  
  if(any(!is.numeric(upper) & !is.na(upper))) {
    stop("upper should consist of numerics and/or NAs.")
  }
  
  if(any(!is.numeric(trunclower))) {
    stop("trunclower should consist of numerics.")
  }
  
  if(any(!is.numeric(truncupper))) {
    stop("truncupper should consist of numerics.")
  }
  
  
  # Check inequalities
  if(!all(is.na(lower)) & !all(is.na(upper))) {
    if(any(lower>upper)) {
      stop("lower should be smaller than (or equal to) upper.")
    }
  }
  
  
  if(any(trunclower>truncupper)) {
    stop("trunclower should be smaller than (or equal to) truncupper.")
  }
  
  
  if(any(!is.finite(trunclower) & !is.finite(truncupper))) {
    stop("trunclower and truncupper cannot be both infinite.")
  }
  
  if(!all(is.na(lower))) {
    if(any(trunclower>lower)) {
      stop("trunclower should be smaller than (or equal to) lower.")
    }
  }
  
  
  if(!all(is.na(upper))) {
    if(any(truncupper<upper)) {
      stop("truncupper should be larger than (or equal to) cupper.")
    }
  }
  
  
  ##
  
  # Check input for eps
  if (!is.numeric(eps) | length(eps)>1) stop("eps should be a numeric of length 1.")
  if (eps<=0) stop("eps should be strictly positive.")
  
  # Check input for beta_tol
  if (!is.numeric(beta_tol) | length(beta_tol)>1) stop("beta_tol should be a numeric of length 1.")
  if (beta_tol>1 | beta_tol<=0) stop("beta_tol should be in (0,1].")
  
  # Check input for maxiter
  if (!is.numeric(maxiter) | length(maxiter)>1) stop("maxiter should be a numeric of length 1.")
  if (maxiter<1) stop("maxiter should be at least 1.")

}



## Tune the initialising parameters M and s using a grid search over the supplied parameter ranges
.MEtune <- function(lower, upper = lower, trunclower = 0, truncupper = Inf, M = 10, s = 1, 
                    nCores = detectCores(), criterium = "AIC", reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
 
  ######
  # Check input
  
  .ME_checkInput(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper,
                 eps=eps, beta_tol=beta_tol, maxiter=maxiter)
 

  tuning_parameters = expand.grid(M, s)
  
  if(nCores==1) {
    
    i <- 1
 
    all_model <- foreach(i = 1:nrow(tuning_parameters), 
                         .export=c(".ME_initial", ".ME_loglikelihood", ".ME_u_z", ".ME_c_z", ".ME_expected_c", ".ME_T", ".theta_nlm", ".ME_em", ".ME_shape_adj", ".ME_shape_red", ".ME_fit"), 
                         .errorhandling = 'pass') %do% {
                           
      suppressWarnings(.ME_fit(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                               M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], 
                               criterium=criterium, reduceM=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter))
    }


  } else {
    
    
    #######
    # Original code
    
    cl <- makePSOCKcluster(nCores)
    registerDoParallel(cl)  

    i <- 1
    all_model <- foreach(i = 1:nrow(tuning_parameters), 
                         .export=c(".ME_initial", ".ME_loglikelihood", ".ME_u_z", ".ME_c_z", ".ME_expected_c", ".ME_T", ".theta_nlm", ".ME_em", ".ME_shape_adj", ".ME_shape_red", ".ME_fit"), 
                         .errorhandling = 'pass') %dopar% {
                           
                           .ME_fit(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper, 
                                   M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], 
                                   criterium=criterium, reduceM=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter)
    }
    stopCluster(cl) 
  }
  
  # Initial values for s and M and obtained values for information criterium and M
  # Return NA when an error occured in fitting
  f1 <- function(x) {
    if (exists(criterium,x)) {
      with(x, get(criterium))
    } else {
      NA
    }
  } 
  crit <- sapply(all_model, f1)
  
  f2 <- function(x) {
    if (exists("M",x)) {
      x$M
    } else {
      NA
    }
  } 
  M_res <-  sapply(all_model, f2)
  performances <- data.frame(tuning_parameters[,1], tuning_parameters[,2], crit, M_res)
  colnames(performances) = c('M_initial', 's', criterium, 'M')
  
  # Select model with lowest IC
  best_index <- which(crit==min(crit,na.rm=TRUE))[1]
  best_model <- all_model[[best_index]]  
  
  list(best_model = best_model, performances = performances, all_model = all_model)
}

## Density function

.ME_density <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, log = FALSE) {
  
  f <- outer(x, shape, dgamma, scale = theta)
  d <- rowSums(t(t(f)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    d <- d / (.ME_cdf(truncupper, theta, shape, alpha) - .ME_cdf(trunclower, theta, shape, alpha)) * ((trunclower <= x) & (x <= truncupper))
  }
  if(log){
    d <- log(d)
  }
  d
}

## Cumulative distribution function

.ME_cdf <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, lower.tail = TRUE, log.p = FALSE) {  
  
  cdf <- outer(x, shape, pgamma, scale=theta)
  p <- rowSums(t(t(cdf)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    l <- .ME_cdf(trunclower, theta, shape, alpha)
    u <- .ME_cdf(truncupper, theta, shape, alpha)
    p <- ((p - l) / (u - l)) ^ {(x <= truncupper)} * (trunclower <= x)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }
  p
} 



# Value-at-Risk (VaR) or quantile function
# p can be a vector
.ME_VaR <- function(p, theta, shape, alpha, trunclower = 0, truncupper = Inf, interval = NULL, start = NULL) {
  
  pl <- length(p)
  VaR <- numeric(pl)
  
  start0 <- start
  interval0 <- interval
  
  # Analytical solution if single shape
  if(length(shape) == 1){
    Ft <- pgamma(trunclower, shape = shape, scale = theta)
    FT <- pgamma(truncupper, shape = shape, scale = theta)
    VaR <- qgamma(p*(FT-Ft)+Ft, shape = shape, scale = theta)
    
  } else {
    
    # Sort p and keep indices
    p_sort <- sort(p, index.return=TRUE)
    p <- p_sort$x
    
    # Move backwards and use previous estimate as starting value and 
    # upper bound for interval when not provided
    for (i in pl:1) {

      if(is.null(interval0)) {
        
        if(trunclower == 0 & truncupper == Inf) {
          interval <- c(qgamma(p[i], shape = min(shape), scale = theta), 
                        qgamma(p[i], shape = max(shape), scale = theta))
          # Possible better upper bound
          if (!is.na(VaR[i+1])) {
            interval[2] <- min(VaR[i+1], interval[2])
          } 
          
        } else { 
          interval <- c(trunclower, min(truncupper, trunclower + qgamma(p[i], shape = max(shape), scale = theta)))
          # Possible better upper bound
          if (!is.na(VaR[i+1])) {
            interval[2] <- min(VaR[i+1], interval[2])
          }
        }
        
      }
      
      if(is.null(start0)) {
        
        if (!is.na(VaR[i+1])) {
          start <- VaR[i+1]
        } else {
          start <- qgamma(p[i], shape = shape[which.max(alpha)], scale = theta)
        }
        
      }
      
      if(p[i]==1){
        # Fixed for truncation case
        VaR[i] <- truncupper
      } else {
        objective <- function(x){return(10000000*(.ME_cdf(x, theta, shape, alpha, trunclower, truncupper)-p[i])^2)}    
        VaR_nlm <- nlm(f = objective, p = start)
        VaR_optimise <- optimise(f = objective, interval = interval)
        VaR[i] <- ifelse(VaR_nlm$minimum < VaR_optimise$objective, VaR_nlm$estimate, VaR_optimise$minimum)    
      
        if(objective(VaR[i])>1e-06){ # in case optimisation fails, retry with more different starting values
          # Fix warnings
          estimate <- NULL
          minimum <- NULL
          
          alpha <- alpha[order(shape)]
          shape <- shape[order(shape)]
          VaR_nlm <-  vector("list", length(shape))
          VaR_optimise <-  vector("list", length(shape))
          # Fixed for truncation case
          interval <- c(trunclower, pmin(truncupper, trunclower + qgamma(p[i], shape, scale = theta)))
          # Possible better upper bound
          if (!is.na(VaR[i+1])) {
            interval[2] <- min(VaR[i+1], interval[2])
          }
          for(j in 1:length(shape)){
            VaR_nlm[[j]] <- nlm(f = objective, p = qgamma(p[i], shape = shape[j], scale = theta))    
            VaR_optimise[[j]] <- optimise(f = objective, interval = interval[c(1, j+1)])
          }
          VaR_nlm <- sapply(VaR_nlm, with, estimate)[which.min(sapply(VaR_nlm, with, minimum))]
          VaR_optimise <- sapply(VaR_optimise, with, minimum)[which.min(sapply(VaR_optimise, with, objective))]
          VaR[i] <- ifelse(objective(VaR_nlm) < objective(VaR_optimise), VaR_nlm, VaR_optimise)  
        }  
      }  

    }

    # Re-order VaR again
    VaR[p_sort$ix] <- VaR
  }
  
  
  return(VaR)  
}


## Value-at-Risk (VaR) or quantile function (old function)
.ME_VaR_old <- function(p, theta, shape, alpha, trunclower = 0, truncupper = Inf, interval = NULL, start = NULL) {
  
  if(is.null(interval)) {
    if(trunclower == 0 & truncupper == Inf){
      interval <- c(qgamma(p, shape = min(shape), scale = theta), 
                    qgamma(p, shape = max(shape), scale = theta))
      
    } else{ 
      interval <- c(trunclower, min(truncupper, trunclower + qgamma(p, shape = max(shape), scale = theta)))
    }
  }
  
  if(is.null(start)) {
    start <- qgamma(p, shape = shape[which.max(alpha)], scale = theta)
  }
  
  if(p==1){
    # Fixed for truncation case
    return(truncupper) 
  }    
  if(length(shape) == 1 & trunclower == 0 & truncupper == Inf){
    VaR <- qgamma(p, shape = shape, scale = theta)
    return(VaR)
  } else{
    objective <- function(x){return(10000000*(.ME_cdf(x, theta, shape, alpha, trunclower, truncupper)-p)^2)}    
    VaR_nlm <- nlm(f = objective, p = start)
    VaR_optimise <- optimise(f = objective, interval = interval)
    VaR <- ifelse(VaR_nlm$minimum < VaR_optimise$objective, VaR_nlm$estimate, VaR_optimise$minimum)    
    
    if(objective(VaR)>1e-06){ # in case optimization fails, retry with more different starting values
      alpha <- alpha[order(shape)]
      shape <- shape[order(shape)]
      VaR_nlm <-  vector("list", length(shape))
      VaR_optimise <-  vector("list", length(shape))
      # Fixed for truncation case
      interval <- c(trunclower, pmin(truncupper, trunclower + qgamma(p, shape, scale = theta)))
      for(i in 1:length(shape)){
        VaR_nlm[[i]] <- nlm(f = objective, p = qgamma(p, shape = shape[i], scale = theta))    
        VaR_optimise[[i]] <- optimise(f = objective, interval = interval[c(1, i+1)])
      }
      VaR_nlm <- sapply(VaR_nlm, with, estimate)[which.min(sapply(VaR_nlm, with, minimum))]
      VaR_optimise <- sapply(VaR_optimise, with, minimum)[which.min(sapply(VaR_optimise, with, objective))]
      VaR <- ifelse(objective(VaR_nlm) < objective(VaR_optimise), VaR_nlm, VaR_optimise)  
    }  
  }
  VaR  
}
# Vectorize function
.ME_VaR_old <- Vectorize(.ME_VaR_old, vectorize.args = c("p", "start"))


## Excess-of-loss reinsurance premium: L xs R

## Excess-of-loss reinsurance premium: C xs R (Retention R, Cover C, Limit L = R+C)

.ME_XL <- function(R, C, theta, shape, alpha) { 
  
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  coeff <- rep(0, M)
  for(n in 1:M){
    coeff[n] <- sum( alphas[0:(M-n)+n] * (0:(M-n)+1) * pgamma(C, shape = 0:(M-n)+2, scale = theta) )
  }
  if(C == Inf){
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) )
  }else{
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) ) +  C *.ME_cdf(R+C, theta, shape, alpha, lower.tail = FALSE)    
  }
  XL  
}
