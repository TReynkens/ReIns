


# Fit splicing of mixed Erlang and (truncated) Pareto
SpliceFitHill <- function(X, const, M = 10, s = 1:10, trunclower = 0,
                          EVTtruncation = FALSE, ncores = NULL) {
  
  n <- length(X)
  
  # Check input for const
  if (any(const<0) | any(const>=1)) {
    stop("const should be a vector of numbers between 0 and 1.")
  }

  if (is.unsorted(const, strictly=TRUE)) {
    stop("const should be a strictly increasing vector.")
  }
  
  l <- length(const)
  
  if (l>1 & any(const==0)) {
    stop("const cannnot have a zero-element.")
  }
  

  if (length(EVTtruncation)!=l & length(EVTtruncation)!=1) {
    stop("EVTtruncation should have length 1 or the same length as const.")
    
  } else if (l!=1 & length(EVTtruncation)==1) {
    EVTtruncation <- rep(EVTtruncation,l)
  }
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  # Determine values of splicing points
    
  Xsort <- sort(X)
  
  # Number of points larger than splicing points
  k_init <- floor((1-const)*n)
  
  # Splicing points
  tvec <- numeric(l)
  tvec <-  Xsort[n-k_init]

  
  # Update const
  const <- 1-k_init/n
  
  # Number of points in splicing part larger than splicing point
  kvec <- k_init
  
  if (l>1) {
    for(i in (l-1):1) {
      kvec[i] <- kvec[i] - sum(kvec[(i+1):l])
    } 
  }
  t1 <- tvec[1]
   
  
  # Mixing Erlang part
  MEind <- (X<=t1) 
  
  # Upper truncated at threshold t
  fit_tune <- MEtune(lower=X[MEind], upper=X[MEind], trunclower=trunclower, truncupper=t1,
                      M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
  MEfit <- fit_tune$best_model
  
  # EVT part
  EVTfit <- list()
  
  type <- character(l)
  
  EVTfit$gamma <- numeric(l)
  EVTfit$endpoint <- rep(Inf, l)

  # Splice parts
  for (i in 1:l) {
    
    # Endpoint for last splicing is Inf
    if (i==l) {
      tt <- Inf
    } else {
      tt <- tvec[i+1]
    }
    
    if (EVTtruncation[i]) {
      res <- trHill(X[X<=tt])
      resDT <- trDT(X[X<=tt], gamma=res$gamma)
      resEndpoint <- trEndpoint(X[X<=tt], gamma=res$gamma, DT=resDT$DT)
      EVTfit$gamma[i] <- res$gamma[res$k==kvec[i]]
      EVTfit$endpoint[i] <- resEndpoint$Tk[resEndpoint$k==kvec[i]]
      type[i] <- "trHill"
      
    } else {
      res <- Hill(X[X<=tt])
      EVTfit$gamma[i] <- res$gamma[res$k==kvec[i]]
      type[i] <- "Hill"
    }
    
  }
  
  return( list(MEfit=MEfit, EVTfit=EVTfit, t=tvec, trunclower=trunclower, const=const, type=type) )
  
}






# Fit splicing of mixed Erlang and Pareto for right or interval censoring
SpliceFitcHill <- function(Z, I = Z, censored, const, M = 10, s = 1:10, trunclower = 0,
                            EVTtruncation = FALSE, ncores = NULL) {
  
  # Check if interval censoring is present
  interval <- !(all(Z==I))
  
  # Check input for const
  if (const<=0 | const>=1) {
    stop("const should be a number strictly between zero and 1.")
  }
  
  if (length(const)>1) {
    stop("const should be a single number.")
  }
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  n <- length(Z)
  k <- floor((1-const)*n)

  Zsort <- sort(Z)
  t <- Zsort[n-k]
  
  if (length(censored)!=n) {
    stop("Z and censored should have the same length.")
  }
  
  if (interval) {
    # Interval censoring
    
    if (length(Z)!=length(I)) {
      stop("Z and I should have the same length.")
    }
    
    # Mixing Erlang part
    MEind <- (Z<=t) 
    
    # Upper truncated at threshold t
    fit_tune <- MEtune(lower=Z[MEind], upper=pmin(I[MEind],t), trunclower=trunclower, truncupper=t,
                       M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
    MEfit <- fit_tune$best_model
    
    # EVT part
    EVTfit <- list()
    if (EVTtruncation) {
      res <- trciHill(Z, I=I, censored=censored)
      resDT <- trciDT(Z, I, censored=censored, gamma1=res$gamma1)
      resEndpoint <- trciEndpoint(Z, I, censored=censored, gamma1=res$gamma1, DT=resDT$DT)
      EVTfit$gamma1 <- res$gamma1[res$k==k]
      EVTfit$endpoint <- resEndpoint$Tk[resEndpoint$k==k]
      type <- "trciHill"
    } else {
      res <- ciHill(Z, I=I, censored=censored)
      EVTfit$gamma1 <- res$gamma1[res$k==k]
      type <- "ciHill"
    }
    
    # const is 1-k/n but use Turnbull for interval censored data
    const <- Turnbull(t,Z,I,censored)
    
  } else {
    # Right censoring
    
    # Mixing Erlang part
    MEind <- (Z<t)
    upper <- Z
    upper[censored] <- t
    
    # Upper truncated at threshold t
    fit_tune <- MEtune(lower=Z[MEind], upper=upper[MEind], trunclower=trunclower, truncupper=t,
                       M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
    MEfit <- fit_tune$best_model
    
    # cHill part
    res <- cHill(Z, censored=censored)
    EVTfit <- list()
    EVTfit$gamma1 <- res$gamma1[res$k==k]
    
    # const is k/n but use Kaplan-Meier for right censored data
    const <- KaplanMeier(t,Z,censored)
    
    type <- "cHill"
  }
  
   return( list(MEfit=MEfit, EVTfit=EVTfit, t=t, trunclower=trunclower, const=const, type=type) )
}



# Fit splicing of mixed Erlang and GPD (POT)
SpliceFitGPD <- function(X, const, M = 10, s = 1:10, trunclower = 0, ncores = NULL) {

  n <- length(X)
  
  # Check input for const
  if (any(const<0) | any(const>=1)) {
    stop("const should be a vector of numbers between 0 and 1.")
  }
  
  if (is.unsorted(const, strictly=TRUE)) {
    stop("const should be a strictly increasing vector.")
  }
  
  l <- length(const)
  
  if (l>1 & any(const==0)) {
    stop("const cannnot have a zero-element.")
  }
  
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  # Determine values of splicing points
  
  Xsort <- sort(X)
  
  # Number of points larger than splicing points
  k_init <- floor((1-const)*n)
  
  # Splicing points
  tvec <- numeric(l)
  tvec <-  Xsort[n-k_init]
  
  
  # Update const
  const <- 1-k_init/n
  
  # Number of points in splicing part larger than splicing point
  kvec <- k_init
  
  if (l>1) {
    for(i in (l-1):1) {
      kvec[i] <- kvec[i] - sum(kvec[(i+1):l])
    } 
  }
  t1 <- tvec[1]
  
  
  # Mixing Erlang part
  MEind <- (X<=t1) 
  
  # Upper truncated at threshold t
  fit_tune <- MEtune(lower=X[MEind], upper=X[MEind], trunclower=trunclower, truncupper=t1,
                     M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
  MEfit <- fit_tune$best_model
  
  # EVT part
  EVTfit <- list()
  
  type <- rep("GPD", l)
  
  EVTfit$gamma <- numeric(l)
  EVTfit$sigma <- numeric(l)
  EVTfit$endpoint <- rep(Inf, l)
  
  
  # Splice parts
  for (i in 1:l) {
    
    # Endpoint for last splicing is Inf
    if (i==l) {
      tt <- Inf
    } else {
      tt <- tvec[i+1]
    }
    
    POTdata <- X[X>tvec[i] & X<=tt]-tvec[i]
    res <- GPDfit(POTdata)
    EVTfit$gamma[i] <- res[1]
    EVTfit$sigma[i] <- res[2]
    
  }

  return( list(MEfit=MEfit, EVTfit=EVTfit, t=tvec, trunclower=trunclower, const=const, type=type))
}

########################################################################

# Splicing PDF
SplicePDF <- function(x, splicefit) {
  
  d <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<=tvec[1])
  
  # Case x<=t
  d[ind] <- const[1] * ME_density(x[ind], shape = MEfit$shape, alpha = MEfit$beta, 
                       theta = MEfit$theta, trunclower = trunclower, truncupper = tvec[1])
  
  # Case x>t
  for (i in 1:l) {
    
    # Next splicing point (Inf for last part)
    tt <- ifelse(i==l, Inf, tvec[i+1])
    
    # Index for all observations in i-th EVTpart
    ind <- x>tvec[i] & x<=tt
    
    # Constant corresponding to next splicing part
    # (1 for last splicing part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Endpoint of splicing part, min of next splicing point and
    # endpoint from Pareto
    e <- min(tt, EVTfit$endpoint[i])

    if (splicefit$type[i]=="GPD") {
      d[ind] <- dtgpd(x[ind], mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e) * (cconst-const[i])
      
    } else if (type[i] %in% c("Hill","cHill","ciHill","trHill","trciHill")) {
      d[ind] <- dtpareto(x[ind], shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e) * (cconst-const[i])
      
    } else {
      stop("Invalid type.")
    }
    # PDF is 0 after endpoint
    d[x>EVTfit$endpoint[i]] <- 0
  }
  
  return(d)
}


# Splicing CDF 
SpliceCDF <- function(x, splicefit) {
  
  p <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<=tvec[1])
  
  # Case x<=t
  p[ind] <- const[1] * ME_cdf(x[ind], shape = MEfit$shape, alpha = MEfit$beta, 
                theta = MEfit$theta, trunclower = trunclower, truncupper = tvec[1])

  # Case x>t
  for (i in 1:l) {
    
    # Next splicing point (Inf for last part)
    tt <- ifelse(i==l, Inf, tvec[i+1])
    
    # Index for all observations in i-th EVTpart
    ind <- x>tvec[i] & x<=tt
    
    # Constant corresponding to next splicing part
    # (1 for last splicing part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Endpoint of splicing part, min of next splicing point and
    # endpoint from Pareto
    e <- min(tt, EVTfit$endpoint[i])
    
    if (splicefit$type[i]=="GPD") {
      # Note that c +F(x)*(1-c) = 1-(1-c)*(1-F(x))
      p[ind] <- const[i] + ptgpd(x[ind], mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e) * (cconst-const[i])
    
    } else if (type[i] %in% c("Hill","cHill","ciHill","trHill","trciHill")) {
      p[ind] <- const[i] + ptpareto(x[ind], shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e) * (cconst-const[i])
    
    } else {
      stop("Invalid type.")
    }
    
    # CDF is 1 after endpoint
    p[x>EVTfit$endpoint[i]] <- 1
  }
  
  return(p)
}



# Splicing quantiles
SpliceQuant <- function(p, splicefit) {
  
  # Check input
  if (!is.numeric(p)) {
    stop("p should be numeric.")
  }
  
  if (any(p<0) | any(p>1)) {
    stop("All elements of p should be in [0,1].")
  }
  
  q <- numeric(length(p))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)

  ind <- (p<const[1])
  
  if (any(ind)) {
    # Quantiles of ME part
    q[ind] <- ME_VaR(p[ind]/const[1], shape = MEfit$shape, alpha = MEfit$beta, 
                     theta = MEfit$theta, trunclower=trunclower, truncupper=tvec[1], 
                     interval=c(trunclower,tvec[1])) 
  }

  
  # Quantiles of EVT part
  
  for (i in 1:l) {
    
    # Next splicing point (Inf for last part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Index for all probabilities in i-th EVTpart
    ind <- p>=const[i] & p<cconst
    
    tt <- ifelse(i==l, Inf, tvec[i+1])
    e <- min(EVTfit$endpoint[i], tt)
    
    if (splicefit$type[i]=="GPD") {
      q[ind] <- qtgpd((p[ind]-const[i])/(cconst-const[i]), mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e)
      
    } else if (splicefit$type[i] %in% c("Hill","cHill","ciHill","trHill","trciHill")) {
      q[ind] <- qtpareto((p[ind]-const[i])/(cconst-const[i]), shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e)
      
    } else {
      stop("Invalid type.")
    }
  }
  
  # Special case
  q[p==1] <- EVTfit$endpoint[l]
  
  return(q)
}


###########################################################################



# Plot of fitted survival function and ECDF estimator + bounds
SpliceECDF <- function(x, X, splicefit, alpha = 0.05, ...) {
  
  plot(x, 1-SpliceCDF(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)

  # ECDF estimator
  fit  <- ecdf(X)
  # Empirical surivival function
  est <- 1-fit(x)
  
  # Confidence bounds using Dvoretzky-Kiefer-Wolfowitz inequality
  # http://stats.stackexchange.com/questions/55500/confidence-intervals-for-empirical-cdf
  n <- length(X)
  eps <- sqrt(1/(2*n)*log(2/alpha))
  lines(x, est, lty=1, col="red")
  lines(x, pmax(est-eps,0), lty=2, col="blue")
  lines(x, pmin(est+eps,1), lty=2, col="blue")
  legend("topright", c("Fitted survival function","Empirical survival function","95% confidence bounds"),
         lty=c(1,1,2), col=c("black","red","blue"))
}

# Plot of fitted survival function and Turnbull estimator + bounds
SpliceTB <- function(x, Z, I = Z, censored, splicefit, alpha = 0.05, ...) {
  
  plot(x, 1-SpliceCDF(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  
  
  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(Z==I, "right", "interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain", conf.int=1-alpha)
  f = stepfun(fit$time,c(1,fit$surv))
  est <- f(x)
  
  lines(x,est,col="red")
  
  # Confidence bounds
  lines(fit$time, fit$lower, col = "blue", lty = 2)
  lines(fit$time, fit$upper, col = "blue", lty = 2)
  lines(x, est, lty=1, col="red")
  legend("topright", c("Fitted survival function","Turnbull estimator","95% confidence bounds"),
         lty=c(1,1,2), col=c("black","red","blue"))
}


# Probability - probability plot with ECDF
SplicePP <- function(x = sort(X), X, splicefit, log = FALSE, ...) {
  
  # ECDF estimator
  fit  <- ecdf(X)
  est <- 1-fit(x)
  
  if (log) {
    ind <- est>0
    plot(-log(est[ind]), -log(1-SpliceCDF(x[ind],splicefit=splicefit)), type="l",
         xlab="-log(Empirical survival probability)",
         ylab="-log(Fitted survival probability)", ...)
  } else {
    plot(est, 1-SpliceCDF(x,splicefit=splicefit), type="l", xlab="Empirical survival probability",
         ylab="Fitted survival probability", ...)
  }
  abline(a=0, b=1)
}

# Probability - probability plot with Turnbull estimator
SplicePP_TB <- function(x = sort(Z), Z, I = Z, censored, splicefit, log = FALSE, ...) {
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  

  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(Z==I, "right", "interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  f = stepfun(fit$time, c(1,fit$surv))
  est <- f(x)
  
  
  if (log) {
    plot(-log(est), -log(1-SpliceCDF(x,splicefit=splicefit)), type="l", xlab="-log(Empirical surv. probs)",
         ylab="-log(Fitted surv. probs)", ...)
  } else {
    plot(est, 1-SpliceCDF(x, splicefit=splicefit), type="l",
         xlab="Empirical survival probability", ylab="Fitted probability", ...)
  }
  abline(a=0,b=1)
}


# Log-log plot with empirical survival function and fitted survival function
SpliceLL <- function(x = sort(X), X, splicefit, ...) {
  
  # ECDF estimator
  fit  <- ecdf(X)

  X <- sort(X)
  plot(log(X), log(1-fit(X)), ylab="log(emp. surv.)", xlab="log(X)", type="p", ...)
  lines(log(x), log(1-SpliceCDF(x, splicefit=splicefit)))
}


# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB <- function(x = sort(Z), Z, I = Z, censored, splicefit, ...) {
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  
  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(Z==I, "right", "interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  f = stepfun(fit$time, c(1,fit$surv))

  Z <- sort(Z)
  plot(log(Z), log(f(Z)), ylab="log(Turnbull surv.)", xlab="log(X)", type="p", ...)
  lines(log(x), log(1-SpliceCDF(x, splicefit=splicefit)))
}

