
# This file contains the procedure to fit the mixed Erlang - Pareto splicing model
# to interval censored data as described in Reynkens et al. (2016) and 
# Section 4.3 in Albrecher et al. (2017).

###########################################################################

# Numerical tolerance
numtol <- .Machine$double.eps^0.5

# Fit splicing model of mixed Erlang and Pareto to interval censored data
.SpliceFiticPareto <- function(L, U, censored, tsplice, M = 3, s = 1:10, trunclower = 0, truncupper = Inf, ncores = NULL, 
                               criterium = c("BIC", "AIC"), reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, cpp = FALSE) {
  
  tsplice <- as.numeric(tsplice)
  
  # Repeat censored if a single number
  if (length(censored) == 1) {
    censored <- rep(censored, length(L))
  }
  
  ######
  # Check input
  .spliceEM_checkInput(L = L, U=U, censored = censored, trunclower = trunclower, tsplice=tsplice, 
                       truncupper = truncupper, eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  
  # Match input argument for criterium
  criterium <- match.arg(criterium)
  
  ##
  # Initial value for pi
  # Faster to already compute it here
  if (requireNamespace("interval", quietly = TRUE)) {
    pi_initial <- 1-.Turnbull2(tsplice, L = L, R = U, censored = censored, trunclower = trunclower, truncupper = truncupper)$surv
    
  } else {
    pi_initial <- 1-Turnbull(tsplice, L = L, R = U, censored = censored, trunclower = trunclower, truncupper = truncupper)$surv
  }
  
  ## Fit model
  fit <- .spliceEM_tune(lower=L, upper=U, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                       M=M, s=s, pi_initial=pi_initial, nCores=ncores, criterium=criterium, reduceM=reduceM, 
                       eps=eps, beta_tol=beta_tol, maxiter=maxiter, cpp=cpp)
  
  # Output as MEfit object
  MEfit <- .MEoutput(fit)
  
  # EVT part
  EVTfit <- EVTfit(gamma=fit$best_model$gamma, endpoint=truncupper, sigma=NULL)
  
  # Return correct type
  type <- c("ME", ifelse(is.finite(truncupper), "TPa", "Pa"))
  
  # Log-likelihood
  loglik <- fit$best_model$loglikelihood
  
  # Compute ICs, parameters: alpha (1 less because sums to 1), r, theta, gamma's, pi's (fit$best_model$pi has length 1 here!)
  df <- 2*MEfit$M-1 + 1 + length(EVTfit$gamma) + length(fit$best_model$pi)
  aic <- -2*loglik + 2*df
  bic <- -2*loglik + log(length(L))*df
  ic <- c(aic, bic)
  names(ic) <- c("AIC", "BIC")
  
  # Return SpliceFit object
  return( SpliceFit(const=fit$best_model$pi, trunclower = trunclower, t=tsplice,  type=type,
                    MEfit=MEfit, EVTfit=EVTfit, loglik=loglik, IC=ic) )
  
}

# Check input arguments for SpliceFiticPareto
.spliceEM_checkInput  <- function(L, U, censored = NULL, trunclower, tsplice, truncupper, eps, beta_tol, maxiter) {
  
  ##
  # Check lengths
  if (length(L) != length(U)) {
    stop("L and U should have equal length.")
  }
  
  
  if (length(trunclower) != 1) {
    stop("trunclower should have length 1.")
  }
  
  if (length(tsplice) != 1) {
    stop("tsplice should have length 1.")
  }
  
  
  if (length(truncupper) != 1 ) {
    stop("truncupper should have length 1.")
  }
  
  ##
  # Check censored
  if (!is.null(censored)) {
    
    if (!is.logical(censored)) {
      stop("censored should be a logical vector.")
    }
    
    if (length(censored) != length(L)) {
      stop("censored should have the same length as L.")
    }
    
    if (sum(censored) == length(L)) {
      stop("Not all observations can be censored.")
    }
    
    if (any(abs(U-L)[!censored] > numtol)) {
      stop("Uncensored observations should have L=U.")
    }
    
    if (any(abs(U-L)[censored] < numtol)) {
      stop("Censored observations cannot have L=U.")
    }
    
  }
  
  
  ##
  # Check data types
  if (any(!is.numeric(L))) {
    stop("L should consist of numerics.")
  }
  
  if (any(!is.numeric(U))) {
    stop("U should consist of numerics.")
  }
  
  if (!is.numeric(trunclower)) {
    stop("trunclower should be numeric.")
  }
  
  if (!is.numeric(tsplice)) {
    stop("tsplice should be numeric.")
  }
  
  if (!is.numeric(truncupper)) {
    stop("truncupper should be numeric.")
  }
  
  
  ##
  # Check inequalities
  if (trunclower < 0) {
    stop("trunclower cannot be strictly negative.")
  }
  
  if (any(L > U)) {
    stop("all elements of L should be smaller than (or equal to) the corresponding elements of U.")
  }
  
  if (trunclower > tsplice) {
    stop("trunclower should be smaller than (or equal to) tsplice.")
  }
  
  if (truncupper <= tsplice) {
    stop("truncupper should be strictly larger than tsplice.")
  }
  
  if (!is.finite(trunclower) ) {
    stop("trunclower has to be finite.")
  }
  
  if (!is.finite(tsplice)) {
    stop("tsplice has to be finite.")
  }
  
  if (any(trunclower > L)) {
    stop("trunclower should be smaller than (or equal to) all elements of L.")
  }
  
  if (any(U > truncupper)) {
    stop("truncupper should be larger than (or equal to) all elements of U.")
  }
  
  ##
  # Check if enough data in all categories
  
  # t^l <= x_i=l_i=u_i <= t <T
  ind1 <- which(!censored & (U <= tsplice))
  # t^l < t < x_i=l_i=u_i <=T
  ind2 <- which(!censored & (U > tsplice))
  # t^l <= l_i < u_i <= t <T
  ind3 <- which(censored & U <= tsplice)
  # t^l < t <= l_i < u_i <=T
  ind4 <- which(censored & L >= tsplice)
  # t^l <= l_i < t < u_i <=T
  #ind5 <- which(censored & L<tsplice & U > tsplice)
  
  if (length(ind1) == 0 & length(ind3) == 0) {
    stop("No data is given for the ME part.")
  }
  
  if (length(ind2) == 0 & length(ind4) == 0) {
    stop("No data is given for the Pareto part.")
  }
  
  ##
  
  # Check input for eps
  if (!is.numeric(eps) | length(eps) > 1) stop("eps should be a numeric of length 1.")
  if (eps <= 0) stop("eps should be strictly positive.")
  
  # Check input for beta_tol
  if (!is.numeric(beta_tol) | length(beta_tol) > 1) stop("beta_tol should be a numeric of length 1.")
  if (beta_tol > 1 | beta_tol <= 0) stop("beta_tol should be in (0,1].")
  
  # Check input for maxiter
  if (!is.numeric(maxiter) | length(maxiter) > 1) stop("maxiter should be a numeric of length 1.")
  if (maxiter < 1) stop("maxiter should be at least 1.")
  
}


## Tune the initialising parameters M and s using a grid search over the supplied parameter ranges
.spliceEM_tune <- function(lower, upper = lower, censored = NULL, trunclower = 0, tsplice, truncupper = Inf, M = 10, s = 1, pi_initial = NULL,
                          nCores = detectCores(), criterium = "BIC", reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, cpp = FALSE) {
  
  # Censoring indicator for each observation
  if (is.null(censored)) {
    censored <- abs(lower-upper) > numtol & !is.na(lower != upper)
  }
  
  # All combinations of M and s
  tuning_parameters <- expand.grid(M, s)
  
  if (nCores == 1) {
    
    # Special case when nCores=1 without any cluster
    
    i <- 1
    
    all_model <- foreach(i = 1:nrow(tuning_parameters), .errorhandling = 'pass') %do% {
      
      suppressWarnings(.spliceEM_fit(lower=lower, upper=upper, censored = censored,
                                     trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                                     M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], pi_initial=pi_initial,
                                     criterium=criterium, reduceM=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, cpp=cpp))
    }
    
    
  } else {
    
    
    #######
    # Original code
    
    cl <- makePSOCKcluster(nCores)
    registerDoParallel(cl)  
    
    i <- 1
    all_model <- foreach(i = 1:nrow(tuning_parameters), .errorhandling = 'pass') %dopar% {
      
      suppressWarnings(.spliceEM_fit(lower=lower, upper=upper, censored = censored, 
                                     trunclower = trunclower, tsplice=tsplice, truncupper = truncupper, 
                                     M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], pi_initial=pi_initial,
                                     criterium=criterium, reduceM=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, cpp=cpp))
    }
    stopCluster(cl) 
  }
  
  
  # Initial values for s and M and obtained values for information criterium and M
  # Return NA when an error occured in fitting
  f1 <- function(x) {
    if (exists(criterium, x)) {
      with(x, get(criterium))
    } else {
      NA
    }
  } 
  crit <- sapply(all_model, f1)
  
  f2 <- function(x) {
    if (exists("M", x)) {
      x$M
    } else {
      NA
    }
  } 

  M_res <-  sapply(all_model, f2)
  performances <- data.frame(tuning_parameters[,1], tuning_parameters[,2], crit, M_res)
  colnames(performances) = c('M_initial', 's', criterium, 'M')
  
  # Select model with lowest IC
  best_index <- which(crit == min(crit, na.rm = TRUE))[1]
  best_model <- all_model[[best_index]]  
  
  return(list(best_model = best_model, performances = performances, all_model = all_model))
}







# Fit splicing model using EM algorithm including shape adjustment and reduction of number of Erlangs
.spliceEM_fit <- function(lower, upper = lower, censored, trunclower = 0, tsplice, truncupper = Inf, M = 10, s = 1, pi_initial = NULL, criterium="AIC", reduceM = TRUE, 
                          eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, cpp = FALSE) {

  # Get initial values
  initial <- .spliceEM_initial(lower=lower, upper=upper, censored = censored, trunclower = trunclower,
                               tsplice=tsplice, truncupper = truncupper, M=M, s=s, pi_initial=pi_initial)
  
  if (cpp) {
    # Use C++ version
    
    # Reduction of the shape combinations
    fit <- .spliceEM_shape_red_cpp(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                                   pi=initial$pi, theta=initial$theta, shape=initial$shape, beta=initial$beta, gamma=initial$gamma, 
                                   criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=FALSE)
    
    # Subsequent adjustment and reduction of the shape combinations
    fit <- .spliceEM_shape_red_cpp(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                                   pi=fit$pi, theta=fit$theta, shape=fit$shape, beta=fit$beta, gamma=fit$gamma,
                                   criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=TRUE)
    
  } else {
    # Use R version
    
    # Reduction of the shape combinations
    fit <- .spliceEM_shape_red(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                               pi=initial$pi, theta=initial$theta, shape=initial$shape, beta=initial$beta, gamma=initial$gamma, 
                               criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=FALSE)
    
    # Subsequent adjustment and reduction of the shape combinations
    fit <- .spliceEM_shape_red(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                               pi=fit$pi, theta=fit$theta, shape=fit$shape, beta=fit$beta, gamma=fit$gamma,
                               criterium=criterium, improve=reduceM, eps=eps, beta_tol=beta_tol, maxiter=maxiter, adj=TRUE)
  }

  # Return fitted parameters of splicing model
  return(list(pi = fit$pi, alpha = fit$alpha, beta = fit$beta, shape = fit$shape, theta = fit$theta, gamma = fit$gamma, 
              loglikelihood = fit$loglikelihood, AIC=fit$AIC, BIC=fit$BIC, M = fit$M, M_initial = M, s = s)) 
}

# Initial values for parameters. 
# Initial value for pi can be provided. When NULL (default), it will be computed using the Turnbull estimator.
.spliceEM_initial <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf, M = 10, s = 1, pi_initial = NULL) {
  
  ##
  # Initial value for pi
  if (is.null(pi_initial)) {
    # Use Turnbull estimator
    
    if (requireNamespace("interval", quietly = TRUE)) {
      pi_initial <- 1-.Turnbull2(tsplice, L=lower, R=upper, censored = censored, trunclower = trunclower, truncupper = truncupper)$surv
      
    } else {
      pi_initial <- 1-Turnbull(tsplice, L=lower, R=upper, censored = censored, trunclower = trunclower, truncupper = truncupper)$surv
    }
    
  } else {
    # Use provided value
    pi <- pi_initial
  }

  ##
  # Initial values for ME parameters, see Section 4.1 in Verbelen et al. (2016)
  
  # Data for initial step: take mean for interval censored data
  uc <- lower[!censored]
  ic <- (lower[censored] + upper[censored]) / 2
  # Right censored
  ind <- which(upper[censored] == truncupper)
  ic[ind] <- (lower[censored])[ind]
  # Left censored
  ind <- which(lower[censored] == trunclower)
  ic[ind] <- (upper[censored])[ind]
  # Combine
  initial_data <- c(uc, ic)
  # Remove initial_data equal to 0 since they cause first shape to be 0
  initial_data[initial_data == 0] <- NA
  
  initial_dataME <- initial_data[initial_data <= tsplice]
  
  # Initial value of theta using spread factor s
  theta <- max(initial_dataME, na.rm = TRUE) / s
  
  # Initial shapes as quantiles
  shape <- unique(ceiling(quantile(initial_dataME, probs = seq(0, 1, length.out = M), na.rm = TRUE) / theta))
  
  # Initial value for alpha's
  alpha <- rep(0, length(shape))
  alpha[1] <- sum(initial_dataME <= shape[1]*theta)
  if (length(shape) > 1) {
    for (i in 2:length(shape)) {
      alpha[i] <- sum(initial_dataME <= shape[i]*theta & initial_dataME > shape[i-1]*theta)
    }
  }
  
  # Keep strictly positive alpha's and corresponding shapes
  shape <- shape[alpha > 0]
  alpha <- alpha[alpha > 0]/sum(alpha)
  
  # alpha to beta
  t_probs <- pgamma(tsplice, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
  beta <- alpha * t_probs / sum(alpha*t_probs)
  
  ##
  # Initial value for gamma
  
  # Hill estimator applied to initial data
  gamma <- .Hillinternal(initial_data, threshold=tsplice)
  
  return(list(pi=pi, gamma = gamma, theta=theta, shape = shape, alpha = alpha, beta=beta))
}


## Reduction of M based on an information criterium: AIC and BIC implemented (4.3 in Verbelen et al. 2015)
.spliceEM_shape_red <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf, pi, theta, shape, beta, gamma, 
                                criterium = "AIC", improve=TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, adj = TRUE) {
  
  # Fit model with initial M
  n <- length(lower)
  
  if (adj) {
    # Adjust shapes
    fit <- .spliceEM_shape_adj(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                               pi=pi, theta=theta, shape = shape, beta=beta, gamma = gamma, 
                               eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  } else {
    # Do not adjust shapes => apply EM directly
    fit <- .spliceEM_fit_raw(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                             pi=pi, theta=theta, shape = shape, beta=beta, gamma = gamma, 
                             eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  }
  
  loglikelihood <- fit$loglikelihood
  IC <- fit[[criterium]]
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta   
  alpha <- fit$alpha
  pi <- fit$pi
  gamma <- fit$gamma
  M <- length(shape)
  
  # Try to improve IC by reducing M
  while(improve && length(shape) > 1) {    
    
    # Remove shape with smallest beta
    min_ind <- which.min(beta)
    new_shape <- shape[-min_ind]
    new_beta <- beta[-min_ind]
    # Ensure that new beta's sum to 1
    new_beta <- new_beta/sum(new_beta)
    
    
    # Fit model with smaller M using new_beta and new_shape as starting values
    if (adj) {
      # Adjust shapes
      fit <- .spliceEM_shape_adj(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                                 pi=pi, theta=theta, shape=new_shape, beta=new_beta, gamma = gamma, 
                                 eps=eps, beta_tol=beta_tol, maxiter=maxiter)
    } else {
      # Do not adjust shapes => apply EM directly
      fit <-  .spliceEM_fit_raw(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                                pi=pi, theta=theta, shape=new_shape, beta=new_beta, gamma = gamma, 
                                eps=eps, beta_tol=beta_tol, maxiter=maxiter)
    }
    
    new_IC <- fit[[criterium]]
    
    # Continue improving if IC is reduced
    if (new_IC < IC) { 
      IC <- new_IC
      loglikelihood <- fit$loglikelihood  
      shape <- fit$shape
      theta <- fit$theta
      beta <- fit$beta  
      alpha <- fit$alpha
      pi <- fit$pi
      gamma <- fit$gamma
      M <- length(shape)
      
    } else {
      improve <- FALSE
    }        
  }
  
  # Compute ICs for final model, 2*length(alpha)-1+1=2*M-1 since we estimate alpha and shape (both have length M),
  # but alpha sums to 1 hence 1 degree of freedom less, and we also estimate theta.
  # Finally, pi and gamma are also estimated.
  df <- 2*length(alpha)-1+1+1+1
  aic <- -2*loglikelihood + 2*df
  bic <- -2*loglikelihood + log(n)*df
  
  return(list(M = M, pi = pi, alpha = alpha, beta = beta, shape = shape, theta = theta, gamma = gamma, 
              loglikelihood = loglikelihood, AIC=aic, BIC=bic)) 
}



# Shape adjustments (4.2 in Verbelen et al. 2015)
.spliceEM_shape_adj <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf, 
                                pi, theta, shape, beta, gamma, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {


  # Fit model with shape equal to the initial value (and M=length(shape))
  fit <- .spliceEM_fit_raw(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice,  truncupper = truncupper, 
                           pi=pi, theta=theta, shape = shape, beta=beta, gamma = gamma, 
                           eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  
  loglikelihood <- fit$loglikelihood
  theta <- fit$theta
  shape <- fit$shape
  beta <- fit$beta
  alpha <- fit$alpha
  pi <- fit$pi
  gamma <- fit$gamma
  M <- length(shape)
  
  # Before and after are the log-likelihoods used in the outer while loop
  before_loglikelihood <- -Inf
  after_loglikelihood <- loglikelihood    
  iteration <- 1
  
  # Improve as long as log-likelihood increases significantly
  while(after_loglikelihood > before_loglikelihood + eps) {   
    
    before_loglikelihood <- after_loglikelihood
    
    # Try changing the shapes (step 1), loop over all shapes
    for(i in M:1) {
      
      improve <- TRUE
      
      # Increase i-th shape as long as improvement and not larger than next shape
      while( improve && (i == M || ifelse(i <= length(shape), shape[i] < shape[i+1]-1, FALSE)) ) {
        
        # Increase i-th shape by 1
        new_shape <- shape
        new_shape[i] <- new_shape[i]+1 
        fit <- .spliceEM_fit_raw(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice,  truncupper = truncupper,
                                 pi=pi, theta=theta, shape=new_shape, beta=beta, gamma = gamma, 
                                 eps=eps, beta_tol=beta_tol, maxiter=maxiter)
        new_loglikelihood <- fit$loglikelihood
        if (new_loglikelihood > loglikelihood + eps) {
          loglikelihood <- new_loglikelihood
          pi <- fit$pi
          shape <- fit$shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          gamma <- fit$gamma
          # number of shape combinations might have changed after EM algorithm if beta_tol > 0
          M <- length(shape)
          i <- min(i, M)
          
        } else {
          improve <- FALSE
        }
      }
    }
    
    # Try decreasing the shapes (step 2), loop over all shapes
    for(i in 1:M) {
      
      improve <- TRUE
      
      # Decrease i-th shape as long as improvement and larger than 1 and not smaller than previous shape
      while( improve && ( (i == 1) || ifelse(i <= length(shape), shape[i] > shape[i-1]+1, FALSE) ) && ifelse(i <= length(shape), shape[i] > 1, FALSE)) {
        
        # Decrease i-th shape by 1
        new_shape <- shape
        new_shape[i] <- new_shape[i]-1
        fit <- .spliceEM_fit_raw(lower=lower, upper=upper, censored = censored, trunclower = trunclower, tsplice=tsplice,  truncupper = truncupper,
                                 pi=pi, theta=theta, shape=new_shape, beta=beta, gamma = gamma, 
                                 eps=eps, beta_tol=beta_tol, maxiter=maxiter)
        new_loglikelihood <- fit$loglikelihood
        if (new_loglikelihood > loglikelihood + eps) {
          loglikelihood <- new_loglikelihood
          pi <- fit$pi
          shape <- fit$shape
          theta <- fit$theta
          beta <- fit$beta 
          alpha <- fit$alpha
          gamma <- fit$gamma
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
  
  # Compute ICs
  df <- 2*length(alpha)-1+1+1+1
  aic <- -2*loglikelihood + 2*df
  bic <- -2*loglikelihood + log(length(lower))*df
  
  return(list(pi = pi, alpha = alpha, beta = beta, shape = shape, theta = theta, gamma = gamma, 
              loglikelihood = loglikelihood, AIC=aic, BIC=bic)) 
}


# Fit splicing model using EM algorithm for given M and shapes
# Additional input: eps = numerical tolerance for EM algorithm (stopping criterion),
# beta_tol = numerical tolerence to discard beta's
# maxiter = maximal number of EM iterations
.spliceEM_fit_raw <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf, 
                              pi, theta, shape, beta, gamma, 
                              eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
 
  n <- length(lower)
  M <- length(shape)
  
  # Separate uncensored and censored observations
  uncensored <- !censored

  # Boolean for having uncensored and censored observations
  #no_censoring <- (sum(uncensored) != 0)  
  #censoring <- (sum(censored) != 0)  
  

  # Indices for the five different cases
  # Some can be numeric(0)!
  
  # t^l <= x_i=l_i=u_i <= t <T
  ind1 <- which(uncensored & (upper <= tsplice))
  # t^l < t < x_i=l_i=u_i <=T
  ind2 <- which(uncensored & (upper > tsplice))
  # t^l <= l_i < u_i <= t <T
  ind3 <- which(censored & upper <= tsplice)
  # t^l < t <= l_i < u_i <=T
  ind4 <- which(censored & lower >= tsplice)
  # t^l <= l_i < t < u_i <=T
  ind5 <- which(censored & lower < tsplice & upper > tsplice)
  
  lower1 <- lower[ind1]
  lower2 <- lower[ind2]
  lower3 <- lower[ind3]; upper3 <- upper[ind3]
  lower4 <- lower[ind4]; upper4 <- upper[ind4]
  lower5 <- lower[ind5]; upper5 <- upper[ind5]
  
  # Truncation probabilities (for ME)
  t_probs <- pgamma(tsplice, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
  
  
  # MLE for alpha, transform beta to alpha
  # alpha_tilde = alpha/(F_ME(t)-F_ME(t^l))
  alpha_tilde <- beta / t_probs
  alpha <- alpha_tilde / sum(alpha_tilde)
  

  # Compute densities and probabilities for observations from the 5 classes
  L <- .spliceEM_densprob(lower1=lower1, lower2=lower2, lower3=lower3, lower4=lower4, lower5=lower5,
                          upper3=upper3, upper4=upper4, upper5=upper5,
                          pi=pi, shape = shape, theta=theta, alpha_tilde=alpha_tilde, gamma = gamma, 
                          trunclower = trunclower, tsplice=tsplice, truncupper = truncupper)
  
  x1_dens_nosum <- L$x1_dens_nosum
  x1_dens <- L$x1_dens
  
  x2_dens <- L$x2_dens

  c3_probs_nosum <- L$c3_probs_nosum
  c3_probs <- L$c3_probs

  c4_probs <- L$c4_probs
  
  c5_probs_nosum <- L$c5_probs_nosum 
  c5_probs <- L$c5_probs

  # Compute log-likelihood
  loglikelihood <- .splice_loglikelihood(x1_dens=x1_dens, x2_dens=x2_dens, c3_probs=c3_probs, c4_probs=c4_probs, c5_probs=c5_probs)
 
  
  iteration <- 1
  old_loglikelihood <- -Inf

  while(loglikelihood - old_loglikelihood > eps & iteration <= maxiter) {
    
    old_loglikelihood <- loglikelihood
    
    
    ############
    # E-step
    ############
    
    if (length(ind5) != 0) {
      
      # Compute P(X_i<=t | t^l<=l_i<=t<u_i, Theta^{(h-1)}) and P(X_i>t | t^l<=l_i<=t<u_i, Theta^{(h-1)})
      probs <- .spliceEM_probs(lower5=lower5, upper5=upper5, trunclower = trunclower, tsplice=tsplice, truncupper = truncupper,
                               pi=pi, theta=theta, shape = shape, alpha = alpha, gamma = gamma)
      
      P1 <- probs$P1
      P2 <- probs$P2
      
    } else {
      P1 <- P2 <- 0
    }

    
    ##
    # pi

    n1_h <- length(ind1) + length(ind3) + sum(P1)
    n2_h <- n - n1_h
    
    ##
    # ME
    
    if (length(ind1) != 0) {
      # Matrix of z_{ij} with i for rows and j for columns
      z1 <- .spliceEM_i_z(x1_dens_nosum=x1_dens_nosum, M=M)
      # sum_{i in S_{i.}} x_i
      E1_ME <- sum(lower1)
        
    } else {
      z1 <- as.matrix(c(0, 0))
	    E1_ME <- 0
    }
    

    if (length(ind3) != 0) {
      z3 <- .spliceEM_iii_z(c3_probs_nosum=c3_probs_nosum, M=M)
      # sum_{i in S_{iii.}} sum_{j=1}^m z3 * E(X_{ij} | Z_{ij}=1, t^l<=l_i<u_i<=t)
      E3_ME <- sum(rowSums(z3 * .spliceEM_Estep_ME_iii(lower3=lower3, upper3=upper3, shape = shape, theta=theta)))
      
    } else {
      z3 <- as.matrix(c(0, 0))
      E3_ME <- 0
    }

    if (length(ind5) != 0) {
      z5 <- .spliceEM_v_z(c5_probs_nosum=c5_probs_nosum, tsplice=tsplice, alpha_tilde=alpha_tilde, 
                          theta=theta, shape = shape, M=M)
      # sum_{i in S_{v.}} P_{1,i} * sum_{j=1}^m z5 * E(X_{ij} | Z_{ij}=1, t^l<=l_i<t<u_i)
      E5_ME <- sum(P1 * rowSums(z5 * .spliceEM_Estep_ME_v(lower5=lower5, tsplice=tsplice, 
                                                          pi=pi, gamma = gamma, shape = shape, theta=theta)))
      
    } else {
      z5 <- as.matrix(c(0, 0))
      E5_ME <- 0
    }

    
    ##
    # gamma
    
    if (length(ind2) != 0) {
      # sum_{i in S_{ii.}} ln(x_i/t)
      E2_Pa <- sum(log(lower2/tsplice))
      
    } else {
      E2_Pa <- 0
    }
    
    
    if (length(ind4) != 0) {
      # sum_{i in S_{iv.}} E( ln(X_i/t) | t^l<t<=l_i<u_i, gamma^{(h-1)})
      E4_Pa <- sum(.spliceEM_Estep_Pa_iv(lower4=lower4, upper4=upper4, gamma = gamma, tsplice=tsplice))
      
    } else {
      E4_Pa <- 0
    }
    
    if (length(ind5) != 0) {
      # sum_{i in S_{v.}} P_{2,i} * E( ln(X_i/t) | t^l<l_i<=t<u_i, gamma^{(h-1)})
      E5_Pa <- sum(P2 * .spliceEM_Estep_Pa_v(upper5=upper5, gamma = gamma, tsplice=tsplice))
      
    } else {
      E5_Pa <- 0
    }

    
    ############
    # M-step
    ############

    ##
    # pi
    pi <- n1_h / n
    
    ##
    # ME
    beta <- (colSums(z1) + colSums(z3) + colSums(z5 * P1)) / n1_h

    # Remove small beta's
    if (min(beta) < beta_tol) {
      shape <- shape[beta > beta_tol]
      beta <- beta[beta > beta_tol]
      beta <- beta / sum(beta)
    }
    
    # Estimate theta
    theta <- exp(nlm(.spliceEM_theta_nlm, log(theta), E1_ME=E1_ME, E3_ME=E3_ME, E5_ME=E5_ME, 
                     n1_h=n1_h, beta=beta, shape = shape, trunclower = trunclower, truncupper=tsplice)$estimate)
 
    ##
    # gamma
    
    gamma_old <- gamma
    
    gamma = 1/n2_h * (E2_Pa + E4_Pa + E5_Pa)
      

    # Solve equation numerically 
    if (is.finite(truncupper)) {
    
      # Minimise to obtain log(gamma)
      .fgamma_trunc <- function(lgammapar, gamma_notrunc) {
        # Ensure that gamma is strictly positive
        gammapar <- exp(lgammapar)
        (gammapar - gamma_notrunc - log(truncupper/tsplice)/((truncupper/tsplice)^(1/gammapar)-1))^2
      }
      
      gamma <- exp(nlm(.fgamma_trunc, p=log(gamma_old), gamma_notrunc=gamma)$estimate)
    } 
    
  
    #######
    iteration <- iteration + 1
    
    # truncation probabilities
    t_probs <- pgamma(tsplice, shape, scale=theta) - pgamma(trunclower, shape, scale=theta)
    
    # MLE for alpha, transform beta to alpha
    alpha_tilde <- beta / t_probs
    alpha <- alpha_tilde / sum(alpha_tilde)
    
    # Compute densities and probabilities for observations from the 5 classes
    L <- .spliceEM_densprob(lower1=lower1, lower2=lower2, lower3=lower3, lower4=lower4, lower5=lower5,
                            upper3=upper3, upper4=upper4, upper5=upper5, 
                            pi=pi, shape = shape, theta=theta, alpha_tilde=alpha_tilde, gamma = gamma, 
                            trunclower = trunclower, tsplice=tsplice, truncupper = truncupper)
    
    x1_dens_nosum <- L$x1_dens_nosum
    x1_dens <- L$x1_dens
    
    x2_dens <- L$x2_dens
    
    c3_probs_nosum <- L$c3_probs_nosum
    c3_probs <- L$c3_probs
    
    c4_probs <- L$c4_probs
    
    c5_probs_nosum <- L$c5_probs_nosum 
    c5_probs <- L$c5_probs
   
    # Compute log-likelihood
    loglikelihood <- .splice_loglikelihood(x1_dens=x1_dens, x2_dens=x2_dens, c3_probs=c3_probs, c4_probs=c4_probs, c5_probs=c5_probs)
  }
  
  # Compute ICs
  df <- 2*length(alpha)-1+1+1+1
  aic <- -2*loglikelihood + 2*df
  bic <- -2*loglikelihood + log(n)*df
  
  return(list(pi = pi, alpha = alpha, beta = beta, shape = shape, theta = theta, gamma = gamma, 
              loglikelihood = loglikelihood,
              iteration = iteration, AIC = aic, BIC = bic))
}


# Compute densities and probabilities for observations from the 5 classes
.spliceEM_densprob <- function(lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5, 
                               pi, shape, theta, alpha_tilde, gamma, 
                               trunclower, tsplice, truncupper) {
  
  n1 <- length(lower1)
  n3 <- length(lower3)
  n5 <- length(lower5)
  
  L <- list()
  
  if (n1 != 0) {
    # Matrix containing Erlang densities for uncensored observations in ME part and all M values of shape
    # Note F_1(x) = sum_{j=1}^M beta_j F_E^t(x) = sum_{j=1}^M beta_j F_E(x)/(F_E(t)-F_E(t^l))
    # and alpha_tilde = beta/(F_E(t)-F_E(t^l))
    L$x1_dens_nosum <- t(t(outer(lower1, shape, dgamma, scale=theta))*alpha_tilde)
    # Combine into density vector (multiply with pi because splicing density)
    L$x1_dens <- pi * rowSums(L$x1_dens_nosum)
    
  } else {
    L$x1_dens_nosum <- as.matrix(c(0, 0))
    L$x1_dens <- numeric(0)
  }

  
  # Densities of uncensored observations in EVT part
  # No problem when ind2=numeric(0)
  L$x2_dens <- (1-pi) * .dtpareto(lower2, gamma = gamma, tsplice=tsplice, truncupper = truncupper)

  if (n3 != 0) {
    # Matrix containing Erlang probabilities for censored observations of type 3 and all M values of shape
    L$c3_probs_nosum <- t(t(outer(upper3, shape, pgamma, scale=theta))*alpha_tilde) - t(t(outer(lower3, shape, pgamma, scale=theta))*alpha_tilde)
    
    # Combine into probability vector (multiply with pi because splicing density): F(u)-F(l)
    L$c3_probs <- pi * rowSums(L$c3_probs_nosum)
    
  } else {
    # Problem when only right-censoring
    L$c3_probs_nosum <- as.matrix(c(0, 0))
    L$c3_probs <- numeric(0)
  }

  
  # Vector of probabilities for uncensored observations of type 4: F(u)-F(l)
  # No problem when ind4=numeric(0)
  L$c4_probs <- pi + (1-pi) * (.ptpareto(upper4, gamma = gamma, tsplice=tsplice, truncupper = truncupper)-
                               .ptpareto(lower4, gamma = gamma, tsplice=tsplice, truncupper = truncupper))


  if (n5 != 0) {
    # Matrix containing Erlang probabilities for censored observations of type 5 and all M values of shape
    L$c5_probs_nosum <- t(t(outer(lower5, shape, pgamma, scale=theta))*alpha_tilde)
    
    # Vector of probabilities for uncensored observations of type 5: F(u)-F(l)
    L$c5_probs <- pi + (1-pi) * .ptpareto(upper5, gamma = gamma, tsplice=tsplice, truncupper = truncupper) - pi * rowSums(L$c5_probs_nosum)
    
  } else {
    # Problem when no intervals over splicing point
    L$c5_probs_nosum <- as.matrix(c(0, 0))
    L$c5_probs <- numeric(0)
  } 

  return(L)
}


  
# Splicing log-likelihood
.splice_loglikelihood <- function(x1_dens, x2_dens, c3_probs, c4_probs, c5_probs) {

  likelihood_contribution <- numeric(0)

  # likelihood contribution (uncensored)
  likelihood_contribution <- c(x1_dens, x2_dens) 

  # likelihood contribution (censored)
  likelihood_contribution <- c(likelihood_contribution, c3_probs, c4_probs, c5_probs) 

  # Compute log-likelihood and put negative contributions equal to -1000
  # ind0 <- which(likelihood_contribution<=0)
  # loglikelihood_contribution <- likelihood_contribution
  # loglikelihood_contribution[!ind0] <- log(likelihood_contribution[!ind0]) 
  # loglikelihood_contribution[ind0] <- -1000
  loglikelihood_contribution <- ifelse(likelihood_contribution > 0, log(likelihood_contribution), -1000)
 
  # log-likelihood
  return(sum(loglikelihood_contribution))
}

####################
# E-step auxiliary functions

# Compute P(X_i<=t | t^l<=l_i<=t<u_i, Theta^{(h-1)}) and P(X_i>t | t^l<=l_i<=t<u_i, Theta^{(h-1)})
# Use pi^{(h-1)}, gamma^{(h-1)}, ...
.spliceEM_probs <- function(lower5, upper5, trunclower, tsplice, truncupper, pi, theta, shape, alpha, gamma) {
  
  # CDF of ME evaluated in l_i (for all i with t^l<=l_i<=t<u_i)
  p1 <- .pME(lower5, trunclower = trunclower, truncupper=tsplice, theta=theta, shape = shape, alpha = alpha)
  
  P1 <- pi*(1-p1) / (pi+(1-pi)*.ptpareto(upper5, gamma = gamma, tsplice=tsplice, truncupper = truncupper) - pi*p1)
  P2 <- 1-P1
  
  return(list(P1=P1, P2=P2))
}


# ^{i.}z_{ij}^{(h)}: posterior probabilities (uncensored)
.spliceEM_i_z <- function(x1_dens_nosum, M) {      

  z1 <- x1_dens_nosum / rowSums(x1_dens_nosum)
  # in case all ^{u}z_{ij}^{(h)} for j=1,...,M are numerically 0
  z1[is.nan(z1)] <- 1/M
  
  return(z1)
}

# ^{iii.}z_{ij}^{(h)}: posterior probabilities (censored)
.spliceEM_iii_z <- function(c3_probs_nosum, M) {      
  
  z3 <- c3_probs_nosum / rowSums(c3_probs_nosum)
  # in case all ^{c}z_{ij}^{(h)} for j=1,...,M are numerically 0
  z3[is.nan(z3)] <- 1/M
  
  return(z3)
}


# ^{v.}z_{ij}^{(h)}: posterior probabilities (censored)
.spliceEM_v_z <- function(c5_probs_nosum, tsplice, alpha_tilde, shape, theta, M) {  

  # Use alpha_tilde since also used in c5_probs_nosum
  pt <- alpha_tilde*pgamma(tsplice, shape = shape, scale=theta)
  z5 <- -sweep(c5_probs_nosum, 2, pt , "-")
  z5 <- z5 / rowSums(z5)
  # in case all ^{cc}z_{ij}^{(h)} for j=1,...,M are numerically 0
  z5[is.nan(z5)] <- 1/M
  
  return(z5)
}


## Expected value of censored observations with l_i<u_i<t for ME
.spliceEM_Estep_ME_iii <- function(lower3, upper3, shape, theta) {  
  
  # E3_ME is a matrix!
  E3_ME <- theta * (outer(upper3, shape+1, pgamma, scale=theta) - outer(lower3, shape+1, pgamma, scale=theta)) / (outer(upper3, shape, pgamma, scale=theta) - outer(lower3, shape, pgamma, scale=theta))
  E3_ME <- t(t(E3_ME)*shape)
  
  # replace numerical 0/0 (NaN) or Inf by correct expected value
  E3_ME <- ifelse(is.nan(E3_ME) | is.infinite(E3_ME), ifelse(outer(lower3, shape*theta, ">"), lower3, upper3), E3_ME)
  
  return(E3_ME)
}

## Expected value of censored observations with l_i<t<u_i for ME
.spliceEM_Estep_ME_v <- function(lower5, tsplice, pi, gamma, shape, theta) {  
  
  # ## E5_ME is a matrix!
  # E5_ME <- .spliceEM_Estep_ME_c(upper3=rep(tsplice,length(lower5)), lower3=lower5, shape = shape, theta=theta) 
  # E5_ME <- sweep(E5_ME, 1, (1-pi) * (upper5^(-1/gamma+1)-tsplice^(-1/gamma+1))/((gamma-1)*tsplice^(-1/gamma)), "+")
  
  E5_ME <- .spliceEM_Estep_ME_iii(lower3=lower5, upper3=rep(tsplice, length(lower5)), shape = shape, theta=theta)
  
  return(E5_ME)
}


# E-step for Pareto: part of E( ln(X_i/t) | t^l<t<=l_i<u_i, gamma^{(h-1)})
.spliceEM_Estep_Pa_iv <- function(lower4, upper4, gamma, tsplice) {  
  
  # Problem if upper4 is infinite
  u4 <- ifelse(is.finite(upper4), upper4/tsplice, 1)
  u4_gamma <- ifelse(is.finite(upper4), (upper4/tsplice)^(-1/gamma), 0)
  
  return( ((log(lower4/tsplice)+gamma)*(lower4/tsplice)^(-1/gamma) - (log(u4)+gamma)*u4_gamma) / ((lower4/tsplice)^(-1/gamma) - u4_gamma) )
}


# E-step for Pareto: part of E( ln(X_i/t) | t^l<l_i<=t<u_i, gamma^{(h-1)})
.spliceEM_Estep_Pa_v <- function(upper5, gamma, tsplice) {  
  
  return( .spliceEM_Estep_Pa_iv(lower4=tsplice, upper4=upper5, gamma = gamma, tsplice=tsplice) )
  
}

####################
# M-step auxiliary functions

# Auxiliary functions used to estimate theta in the M-step

# This function returns log of optimal theta!
.spliceEM_theta_nlm <- function(ltheta, E1_ME, E3_ME, E5_ME, n1_h, beta, shape, trunclower, truncupper) {
  
  theta <- exp(ltheta)
  TT <- .spliceEM_T(trunclower = trunclower, truncupper = truncupper, shape = shape, theta=theta, beta=beta)
  
  return( (theta - ((E1_ME + E3_ME + E5_ME)/n1_h-TT) / sum(beta*shape)) ^ 2 )
}

.spliceEM_T <- function(trunclower, truncupper, shape, theta, beta) {
  
  # Avoid NaN
  if (truncupper == Inf) { 
    
    # Take log first for numerical stability (avoid Inf / Inf)
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(1 - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)   
    
  } else {    
    
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc_log_2 <- shape*log(truncupper)-truncupper/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)-exp(deriv_trunc_log_2)
  }
  
  return(sum(beta*deriv_trunc))
}



#####################
# Fast CDFs and PDFs


# ME CDF (faster version based on .ME_cdf)
.pME <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf) { 
  
  # Matrix with pgamma(x, shape, scale=theta) for elements of x (rows) and elements of shape (columns)
  cdf <- outer(x, shape, pgamma, scale=theta)
  # Compute CDF as product with alpha and then sum per column
  p <- colSums(t(cdf)*alpha)
  
  if (!(trunclower == 0 & truncupper == Inf)) {
    l <- .pME(trunclower, theta=theta, shape = shape, alpha = alpha)
    u <- .pME(truncupper, theta=theta, shape = shape, alpha = alpha)
    p <- ((p - l) / (u - l)) ^ {(x <= truncupper)} * (trunclower <= x)
  }
  
  return(p)
} 

# Pareto PDF
.dpareto <- function(x, gamma, tsplice) {
  
  ifelse(x <= tsplice, 0, 1/(gamma*tsplice) * (x/tsplice)^(-1/gamma-1))
}


# Pareto CDF
.ppareto <- function(x, gamma, tsplice) {
  
  ifelse(x <= tsplice, 0, 1-(x/tsplice)^(-1/gamma))
}

# Truncated Pareto PDF
.dtpareto <- function(x, gamma, tsplice, truncupper) {
  
  ifelse(x <= truncupper, .dpareto(x=x, gamma = gamma, tsplice=tsplice)/.ppareto(x=truncupper, gamma = gamma, tsplice=tsplice), 0)
}

# Truncated Pareto CDF
.ptpareto <- function(x, gamma, tsplice, truncupper) {
  
  ifelse(x <= truncupper, .ppareto(x=x, gamma = gamma, tsplice=tsplice)/.ppareto(x=truncupper, gamma = gamma, tsplice=tsplice), 1)
}
