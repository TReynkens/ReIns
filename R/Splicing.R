

# Check input for const
constCheck <- function(const) {
 
  if (!is.numeric(const)) stop("const should be numeric.")
  
  if (any(const<=0) | any(const>=1)) {
    stop("const should be a vector of numbers in (0,1).")
  }
  
  if (is.unsorted(const, strictly=TRUE)) {
    stop("const should be a strictly increasing vector.")
  }
  
}


#########################################################
# S3 classes


# Check if x is an integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
} 


# MEfit class
MEfit <- function(p, shape, theta, M, M_initial = NULL) {
  
  # Check if arguments are numeric
  if (!is.numeric(p)) stop("p must be numeric.")
  if (!is.numeric(shape)) stop("shape must be numeric.")
  if (!is.numeric(theta)) stop("theta must be numeric.")
  if (!is.numeric(M)) stop("M must be numeric.")
  
  # Check length of M and if it is an integer
  if (length(M)!=1) stop("M must have length 1.")
  if (!is.wholenumber(M) | M<=0) stop("M must be a strictly positive integer.")
  
  # Check length of arguments
  if (length(p)!=M) stop("p must have length M.")
  if (length(shape)!=M) stop("shape must have length M.")
  if (length(theta)!=1) stop("theta must have length 1")
  
  # Check if p is a vector of probabilities that sums to one.
  if (any(p<0)) stop("p must contain positive numbers.")
  if (abs(sum(p)-1)>10^(-16)) stop("p must have sum 1.")
  
  # Make first list
  L <- list(p=p, shape=shape, theta=theta, M=M)
  
  # Optional argument M_initial
  if (!is.null(M_initial)) {
    
    if (!is.numeric(M_initial)) stop("M_initial must be numeric.")
    
    if (length(M_initial)!=1) stop("M_initial must have length 1.")
    
    if (!is.wholenumber(M_initial) | M_initial<=0) {
      stop("M_initial must be a strictly positive integer.")
    }
    
    L$M_initial <- M_initial
  }
  
  # Return final structure of class MEfit
  structure(L, class = "MEfit")
}


# EVTfit class
EVTfit <- function(gamma, endpoint = NULL, sigma = NULL) {
  
  # Default value
  if (is.null(endpoint)) endpoint <- rep(Inf, length(gamma))
  
  # Check if arguments are numeric
  if (!is.numeric(gamma)) stop("gamma must be numeric.")
  if (!is.numeric(endpoint)) stop("endpoint must be numeric.")
  
  
  # Check length of gamma and endpoint
  if (length(gamma)!=length(endpoint)) stop("gamma and endpoint should have equal length.")
  
  # Make first list
  L <- list(gamma=gamma)
  
  # Optional argument sigma
  if (!is.null(sigma)) {
    
    if (!is.numeric(sigma)) stop("sigma must be numeric.")
    
    if (length(sigma)!=length(gamma)) stop("gamma and sigma should have equal length.")
    
    if (any(sigma<=0)) {
      stop("sigma should be strictly positive.")
    }
    
    L$sigma <- sigma
  }
  
  # Add endpoint to list
  L$endpoint <- endpoint
  
  # Return final structure of class EVTfit
  structure(L, class = "EVTfit")
}



# SpliceFit class
SpliceFit <- function(const, trunclower, t, type, MEfit, EVTfit) {
  
  # Check input for const
  constCheck(const)
  
  l <- length(const)
  
  
  # Compute pi
  pi <- c(const[1], diff(const), 1-const[l])

  
  # Check t and trunclower
  if (!is.numeric(trunclower)) stop("trunclower should be numeric.")
  if (length(trunclower)!=1) stop("trunclower should have length 1.")
  
  if(length(t)!=l) stop("const and t should have equal length.")
  
  if (is.unsorted(t, strictly=TRUE)) {
    stop("t should be a strictly increasing vector.")
  }
  
  if (t[1]<=trunclower) stop("trunclower should be strictly smaller than the elements of t.")
  
  
  # Check type
  if (length(type)!=l+1) stop("type should have one element more than const.")
  
  if (!type[1]=="ME") stop("The first argument of type is \"ME\".")
  if (!all(type[-1] %in% c("Pa", "GPD", "tPa", "cPa", "ciPa", "tciPa"))) stop("Invalid type.")
  
  
  # Check MEfit and EVTfit
  if (class(MEfit)!="MEfit") stop("MEfit should be of class MEfit.")
  if (class(EVTfit)!="EVTfit") stop("EVTfit should be of class EVTfit.")
  
  if (length(EVTfit$gamma)!=l) stop("gamma should have the same length as const.")
  
  
  # Make list
  L <- list(const=const, pi=pi, trunclower=trunclower, t=t, type=type, MEfit=MEfit, EVTfit=EVTfit)
  
  # Return final structure of class Splicefit
  structure(L, class = "SpliceFit")
}


# Only keep necessary parts from MEtune output
MEoutput <- function(fit_tune) {
  
  # Obtain best model
  MEfit_old <- fit_tune$best_model
  
  # Make MEfit object
  MEfit <- MEfit(p=MEfit_old$beta, shape=MEfit_old$shape, theta=MEfit_old$theta, 
                 M=MEfit_old$M, M_initial=MEfit_old$M_initial)

  return(MEfit)
}


###############################################################################

# Fit splicing of mixed Erlang and (truncated) Pareto
SpliceFitHill <- function(X, const, M = 3, s = 1:10, trunclower = 0,
                          EVTtruncation = FALSE, ncores = NULL) {
 
  # Check if X is numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  n <- length(X)
  
  # Check input for const
  constCheck(const)
  
  l <- length(const)
  

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

  # Problem when first splicing point smaller than trunclower
  if (trunclower>=tvec[1]) {
    stop("trunclower should be strictly smaller than the first splicing point.")
  }
  
  if (trunclower>min(X)) {
    stop("trunclower should be strictly smaller than all data points.")
  }
  
  
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
  # Output as MEfit object
  MEfit <- MEoutput(fit_tune)
  
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
      type[i] <- "tPa"
      
    } else {
      res <- Hill(X[X<=tt])
      EVTfit$gamma[i] <- res$gamma[res$k==kvec[i]]
      type[i] <- "Pa"
    }
    
  }
  
  # Convert to object of class EVTfit
  EVTfit <- structure(EVTfit, class="EVTfit")
  
  # Return SpliceFit object
  return( SpliceFit(const=const, trunclower=trunclower, t=tvec,  type=c("ME",type), MEfit=MEfit, EVTfit=EVTfit) )
}






# Fit splicing of mixed Erlang and Pareto for right or interval censoring
SpliceFitcHill <- function(Z, I = Z, censored, const, M = 3, s = 1:10, trunclower = 0,
                            EVTtruncation = FALSE, ncores = NULL) {
  
  # Check if Z and I are numeric
  if (!is.numeric(Z)) stop("Z should be a numeric vector.")
  if (!is.numeric(I)) stop("I should be a numeric vector.")
  
  
  # Check if interval censoring is present
  interval <- !(all(Z==I))
  
  # Check input for const
  constCheck(const)
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  n <- length(Z)
  k <- floor((1-const)*n)

  Zsort <- sort(Z)
  t <- Zsort[n-k]
  
  # Problem when splicing point smaller than trunclower
  if (trunclower>=t) {
    stop("trunclower should be strictly smaller than the splicing point.")
  }
  
  if (trunclower>min(Z)) {
    stop("trunclower should be strictly smaller than all data points.")
  }
  
  
  
  
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
    # Output as MEfit object
    MEfit <- MEoutput(fit_tune)

    # EVT part
    EVTfit <- list()
    if (EVTtruncation) {
      res <- trciHill(Z, I=I, censored=censored)
      resDT <- trciDT(Z, I, censored=censored, gamma1=res$gamma1)
      resEndpoint <- trciEndpoint(Z, I, censored=censored, gamma1=res$gamma1, DT=resDT$DT)
      EVTfit$gamma <- res$gamma1[res$k==k]
      EVTfit$endpoint <- resEndpoint$Tk[resEndpoint$k==k]
      type <- "tciPa"
    } else {
      res <- ciHill(Z, I=I, censored=censored)
      EVTfit$gamma <- res$gamma1[res$k==k]
      EVTfit$endpoint <- Inf
      type <- "ciPa"
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
    # Output as MEfit object
    MEfit <- MEoutput(fit_tune)
    
    # cHill part
    res <- cHill(Z, censored=censored)
    EVTfit <- list()
    EVTfit$gamma <- res$gamma1[res$k==k]
    EVTfit$endpoint <- Inf
    
    # const is k/n but use Kaplan-Meier for right censored data
    const <- KaplanMeier(t,Z,censored)
    
    type <- "cPa"
  }
  
  # Convert to object of class EVTfit
  EVTfit <- structure(EVTfit, class="EVTfit")
  
  # Return SpliceFit object
  return( SpliceFit(const=const, trunclower=trunclower, t=t,  type=c("ME",type), MEfit=MEfit, EVTfit=EVTfit) )
}



# Fit splicing of mixed Erlang and GPD (POT)
SpliceFitGPD <- function(X, const, M = 3, s = 1:10, trunclower = 0, ncores = NULL) {

  # Check if X is numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  n <- length(X)
  
  
  # Check input for const
  constCheck(const)
  
  l <- length(const)
  
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
  
  # Problem when first splicing point smaller than trunclower
  if (trunclower>=tvec[1]) {
    stop("trunclower should be strictly smaller than the first splicing point.")
  }
  
  if (trunclower>min(X)) {
    stop("trunclower should be strictly smaller than all data points.")
  }
  
  
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
  # Output as MEfit object
  MEfit <- MEoutput(fit_tune)
  
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

  # Convert to object of class EVTfit
  EVTfit <- structure(EVTfit, class="EVTfit")
  
  # Return SpliceFit object
  return( SpliceFit(const=const, trunclower=trunclower, t=tvec,  type=c("ME",type), MEfit=MEfit, EVTfit=EVTfit) )
}



########################################################################

# Splicing PDF
dSplice <- function(x, splicefit, log = FALSE) {
  
  d <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<tvec[1])
  
  # Case x<=t
  d[ind] <- const[1] * ME_density(x[ind], shape = MEfit$shape, alpha = MEfit$p, 
                       theta = MEfit$theta, trunclower = trunclower, truncupper = tvec[1])
  
  # Case x>t
  for (i in 1:l) {
    
    # Next splicing point (Inf for last part)
    tt <- ifelse(i==l, Inf, tvec[i+1])
    
    # Index for all observations in i-th EVTpart
    ind <- x>=tvec[i] & x<tt
    
    # Constant corresponding to next splicing part
    # (1 for last splicing part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Endpoint of splicing part, min of next splicing point and
    # endpoint from Pareto
    e <- min(tt, EVTfit$endpoint[i])

    # i+1 since type[1]="ME"
    if (splicefit$type[i+1]=="GPD") {
      d[ind] <- dtgpd(x[ind], mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e) * (cconst-const[i])
      
    } else if (type[i+1] %in% c("Pa","cPa","ciPa","tPa","tciPa")) {
      d[ind] <- dtpareto(x[ind], shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e) * (cconst-const[i])
      
    } else {
      stop("Invalid type.")
    }
    # PDF is 0 after endpoint
    d[x>=EVTfit$endpoint[i]] <- 0
  }
  
  if (log) d <- log(d)
  
  return(d)
}


# Splicing CDF 
pSplice <- function(x, splicefit, lower.tail = TRUE, log.p = FALSE) {
  
  p <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<tvec[1])
  
  # Case x<t
  p[ind] <- const[1] * ME_cdf(x[ind], shape = MEfit$shape, alpha = MEfit$p, 
                theta = MEfit$theta, trunclower = trunclower, truncupper = tvec[1])

  # Case x>t
  for (i in 1:l) {
    
    # Next splicing point (Inf for last part)
    tt <- ifelse(i==l, Inf, tvec[i+1])
    
    # Index for all observations in i-th EVTpart
    ind <- x>=tvec[i] & x<tt
    
    # Constant corresponding to next splicing part
    # (1 for last splicing part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Endpoint of splicing part, min of next splicing point and
    # endpoint from Pareto
    e <- min(tt, EVTfit$endpoint[i])
    
    # i+1 since type[1]="ME"
    if (splicefit$type[i+1]=="GPD") {
      # Note that c +F(x)*(1-c) = 1-(1-c)*(1-F(x))
      p[ind] <- const[i] + ptgpd(x[ind], mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e) * (cconst-const[i])
    
    } else if (type[i+1] %in% c("Pa","cPa","ciPa","tPa","tciPa")) {
      p[ind] <- const[i] + ptpareto(x[ind], shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e) * (cconst-const[i])
    
    } else {
      stop("Invalid type.")
    }
    
    # CDF is 1 after endpoint
    p[x>=EVTfit$endpoint[i]] <- 1
  }
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  
  return(p)
}



# Splicing quantiles
qSplice <- function(p, splicefit, lower.tail = TRUE, log.p = FALSE) {
  
  # Check input
  if (!is.numeric(p)) {
    stop("p should be numeric.")
  }
  
  if (log.p) p <- exp(p)
  
  
  if (any(p<0) | any(p>1)) {
    stop("All elements of p should be in [0,1].")
  }
  
  if (!lower.tail) p <- 1-p
  
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
    q[ind] <- ME_VaR(p[ind]/const[1], shape = MEfit$shape, alpha = MEfit$p, 
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
    
    # i+1 since type[1]="ME"
    if (splicefit$type[i+1]=="GPD") {
      q[ind] <- qtgpd((p[ind]-const[i])/(cconst-const[i]), mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e)
      
    } else if (splicefit$type[i+1] %in% c("Pa","cPa","ciPa","tPa","tciPa")) {
      q[ind] <- qtpareto((p[ind]-const[i])/(cconst-const[i]), shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e)
      
    } else {
      stop("Invalid type.")
    }
  }
  
  # Special case
  q[p==1] <- EVTfit$endpoint[l]
  
  return(q)
}

# Random numbers from splicing distribution
rSplice <- function(n, splicefit) {
  
  return(qSplice(runif(n), splicefit=splicefit))
}

###########################################################################



# Plot of fitted survival function and ECDF estimator + bounds
SpliceECDF <- function(x, X, splicefit, alpha = 0.05, ...) {
  
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
  lines(x, pmax(est-eps,0), lty=2, col="blue")
  lines(x, pmin(est+eps,1), lty=2, col="blue")
  legend("topright", c("Fitted survival function","Empirical survival function","95% confidence bounds"),
         lty=c(1,1,2), col=c("black","red","blue"))
}

# Plot of fitted survival function and Turnbull estimator + bounds
SpliceTB <- function(x, Z, I = Z, censored, splicefit, alpha = 0.05, ...) {
  
  plot(x, 1-pSplice(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  
  
  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(all.equal(Z,I), "right", "interval")

  if (type=="interval") {
    fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain", conf.int=1-alpha)
  } else {
    fit  <- survfit(Surv(time=L, event=event, type=type) ~1, conf.type="plain", conf.int=1-alpha)
  }
  
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
    plot(-log(est[ind]), -log(1-pSplice(x[ind],splicefit=splicefit)), type="l",
         xlab="-log(Empirical survival probability)",
         ylab="-log(Fitted survival probability)", ...)
  } else {
    plot(est, 1-pSplice(x,splicefit=splicefit), type="l", xlab="Empirical survival probability",
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
  type <- ifelse(all.equal(Z,I), "right", "interval")
  
  if (type=="interval") {
    fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  } else {
    fit  <- survfit(Surv(time=L, event=event, type=type) ~1, conf.type="plain")
  }
  
  f = stepfun(fit$time, c(1,fit$surv))
  est <- f(x)
  
  
  if (log) {
    plot(-log(est), -log(1-pSplice(x,splicefit=splicefit)), type="l", xlab="-log(Empirical surv. probs)",
         ylab="-log(Fitted surv. probs)", ...)
  } else {
    plot(est, 1-pSplice(x, splicefit=splicefit), type="l",
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
  lines(log(x), log(1-pSplice(x, splicefit=splicefit)))
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
  type <- ifelse(all.equal(Z,I), "right", "interval")
  
  if (type=="interval") {
    fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  } else {
    fit  <- survfit(Surv(time=L, event=event, type=type) ~1, conf.type="plain")
  }
  f = stepfun(fit$time, c(1,fit$surv))

  Z <- sort(Z)
  plot(log(Z), log(f(Z)), ylab="log(Turnbull surv.)", xlab="log(X)", type="p", ...)
  lines(log(x), log(1-pSplice(x, splicefit=splicefit)))
}

