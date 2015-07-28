


# Fit splicing of mixed Erlang and (truncated) Pareto
# If const is non-zero, an extra splicing with a Pareto is included
SpliceFitHill <- function(X, const, const2 = 0, M=20, s=1:10, trunclower = 0,
                          EVTtruncation = FALSE, ncores = NULL) {
  
  n <- length(X)
  
  # Check input for const and const2
  if(const<=0 | const>=1) {
    stop("const should be a number strictly between zero and 1.")
  }
  
  if(const2<0 | const2>=1) {
    stop("const2 should be a number between zero and 1.")
  }
  
  if(const2!=0 & const2<=const) {
    stop("const2 should be 0 or strictly larger than const.")
  }
  
  if(const2>0) {
    # Double splicing
    k2 <- floor((1-const2)*n)
    k1 <- floor((1-const)*n)-k2
  } else {
    # Single splicing
    k1 <- floor((1-const)*n)
  }

  
  # Check input for ncores
  if(is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if(is.na(ncores)) ncores <- 1
  
  Xsort <- sort(X)
  t1 <- Xsort[length(X)-k1]
  if(const2>0) {
    t2 <- Xsort[length(X)-k2]
  } else {
    t2 <- Inf
  }
  
  
  # Mixing Erlang part
  MEind <- (X<=t1) 
  
  # Upper truncated at threshold t
  fit_tune <- MEtune(lower=X[MEind], upper=X[MEind], trunclower=trunclower, truncupper=t1,
                      M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
  MEfit <- fit_tune$best_model
  
  # EVT part
  EVTfit <- list()
  if (EVTtruncation) {
    res <- trHill(X[X<=t2])
    resDT <- trDT(X[X<=t2], gamma=res$gamma)
    resEndpoint <- trEndpoint(X[X<=t2], gamma=res$gamma, DT=resDT$DT)
    EVTfit$gamma <- res$gamma[res$k==k1]
    EVTfit$endpoint <- resEndpoint$Tk[resEndpoint$k==k1]
    type <- "trHill"
  } else {
    res <- Hill(X[X<=t2])
    EVTfit$gamma <- res$gamma[res$k==k1]
    type <- "Hill"
  }
  
  const <- 1-k1/n

  if(const2==0) {
    return(list(MEfit=MEfit,EVTfit=EVTfit,t=t1,trunclower=trunclower,const=const,type=type))
    
  } else {
    # Second splicing part
    res <- Hill(X)
    EVTfit$gamma2 <- res$gamma[res$k==k2]
    
    const2 <- 1-k2/n
    return(list(MEfit=MEfit,EVTfit=EVTfit,t=t1,t2=t2,
                trunclower=trunclower,const=const,const2=const2,type=type))
  }
  
}






# Fit splicing of mixed Erlang and Pareto for right or interval censoring
SpliceFitcHill <- function(Z, I = Z, censored, const, M=20, s=1:10, trunclower = 0,
                            EVTtruncation = FALSE, ncores = NULL) {
  
  # Check if interval censoring is present
  interval <- !(all(Z==I))
  
  # Check input for const
  if(const<=0 | const>=1) {
    stop("const should be a number strictly between zero and 1.")
  }
  
  
  # Check input for ncores
  if(is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if(is.na(ncores)) ncores <- 1
  
  
  n <- length(Z)
  k <- floor((1-const)*n)

  Zsort <- sort(Z)
  t <- Zsort[n-k]
  
  if(length(censored)!=n) {
    stop("Z and censored should have the same length.")
  }
  
  if(interval) {
    # Interval censoring
    
    if(length(Z)!=length(I)) {
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
  
   return(list(MEfit=MEfit,EVTfit=EVTfit,t=t,trunclower=trunclower,const=const,type=type))
}



# Fit splicing of mixed Erlang and GPD (POT)
SpliceFitGPD <- function(X, const, M=20, s=1:10, trunclower = 0, ncores = NULL) {

  # Check input for ncores
  if(is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if(is.na(ncores)) ncores <- 1
  
  # Check input for const
  if(const<=0 | const>=1) {
    stop("const should be a number strictly between zero and 1.")
  }
  
  n <- length(X)
  k <- floor((1-const)*n)
  
  Xsort <- sort(X)
  t <- Xsort[length(X)-k]
  
  # Mixing Erlang part
  MEdata <- X[X<=t]
  #Upper truncated at threshold t
  fit_tune <- MEtune(lower=MEdata, upper=MEdata, trunclower=trunclower, truncupper=t,
                      M=M, s=s, nCores = ncores, criterium="AIC", eps=1e-03, print=FALSE)
  MEfit <- fit_tune$best_model
  
  # POT part
  
  POTdata <- X[X>t]-t
  res <- GPDfit(POTdata)
  POTfit <- list()
  POTfit$gamma <- res[1]
  POTfit$sigma <- res[2]
  
  return(list(MEfit=MEfit,EVTfit=POTfit,t=t,trunclower=trunclower,const=1-k/n,type="GPD"))
}

########################################################################

# Splicing PDF
SplicePDF <- function(x, splicefit) {
  
  d <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  t <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  
  ind <- x<=t
  
  # Case x<=t
  d[ind] <- const * ME_density(x[ind], shape = MEfit$shape, alpha = MEfit$beta, 
                       theta = MEfit$theta, trunclower = trunclower, truncupper = t)
  
  # Case x>t
  if(exists("const2",where=splicefit)) {
    # Double splicing
    
    const1 <- const
    t1 <- t
    t2 <- splicefit$t2
    const2 <- splicefit$const2
    ind2 <- x>t1 & x<=t2
    
    # Case x>t
    if(splicefit$type=="GPD") {
      # Note that c +F(x)*(1-c) = 1-(1-c)*(1-F(x))
      d[ind2] <- dgpd(x[ind2],gamma=EVTfit$gamma,mu=t1,sigma=EVTfit$sigma) * (const2-const1)
    } else if(splicefit$type %in% c("Hill","cHill","ciHill")) {
      d[ind2] <- dpareto(x[ind2],shape=1/EVTfit$gamma,scale=t1) * (const2-const1)
    } else if(splicefit$type %in% c("trHill","trciHill")) {
      d[ind2] <- dtpareto(x[ind2],shape=1/EVTfit$gamma,scale=t1,endpoint=EVTfit$endpoint) * (const2-const1)
      # CDF is 1 after endpoint
      d[x>EVTfit$endpoint] <- 1
    }
    
    ind3 <- x>t2
    d[ind3] <- dpareto(x[ind3],shape=1/EVTfit$gamma2,scale=t2) * (1-const2)
    
  } else {
    
    if(splicefit$type=="GPD") {
      d[!ind] <- (1-const) * dgpd(x[!ind],gamma=EVTfit$gamma,mu=t,sigma=EVTfit$sigma)
    } else if(splicefit$type %in% c("Hill","cHill","ciHill")) {
      d[!ind] <- (1-const) * dpareto(x[!ind], shape=1/EVTfit$gamma, scale=t)
    } else if(splicefit$type %in% c("trHill","trciHill")) {
      d[!ind] <- (1-const) * dtpareto(x[!ind], shape=1/EVTfit$gamma, scale=t, endpoint=EVTfit$endpoint)
      # density is 0 after endpoint
      d[x>EVTfit$endpoint] <- 0
    }
  }

  return(d)
}


# Splicing CDF 
SpliceCDF <- function(x, splicefit) {
  
  p <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  t <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  
  ind <- (x<=t)
  
  # Case x<=t
  p[ind] <- const * ME_cdf(x[ind], shape = MEfit$shape, alpha = MEfit$beta, 
                theta = MEfit$theta, trunclower = trunclower, truncupper = t)

  # Case x>t
  if(exists("const2",where=splicefit)) {
    # Double splicing

    const1 <- const
    t1 <- t
    t2 <- splicefit$t2
    const2 <- splicefit$const2
    ind2 <- x>t1 & x<=t2
    
    # Case x>t
    if(splicefit$type=="GPD") {
      # Note that c +F(x)*(1-c) = 1-(1-c)*(1-F(x))
      p[ind2] <- const1 + pgpd(x[ind2],gamma=EVTfit$gamma,mu=t1,sigma=EVTfit$sigma) * (const2-const1)
    } else if(splicefit$type %in% c("Hill","cHill","ciHill")) {
      p[ind2] <- const1 + ppareto(x[ind2],shape=1/EVTfit$gamma,scale=t1) * (const2-const1)
    } else if(splicefit$type %in% c("trHill","trciHill")) {
      p[ind2] <- const1 + ptpareto(x[ind2],shape=1/EVTfit$gamma,scale=t1,endpoint=EVTfit$endpoint) * (const2-const1)
      # CDF is 1 after endpoint
      p[x>EVTfit$endpoint] <- 1
    }
    
    ind3 <- x>t2
    p[ind3] <- const2 + ppareto(x[ind3],shape=1/EVTfit$gamma2,scale=t2) * (1-const2)
    
    
  } else {
    
    if(splicefit$type=="GPD") {
      # Note that c +F(x)*(1-c) = 1-(1-c)*(1-F(x))
      p[!ind] <- const + pgpd(x[!ind],gamma=EVTfit$gamma,mu=t,sigma=EVTfit$sigma) * (1-const)
    } else if(splicefit$type %in% c("Hill","cHill","ciHill")) {
      p[!ind] <- const + ppareto(x[!ind],shape=1/EVTfit$gamma,scale=t) * (1-const)
    } else if(splicefit$type %in% c("trHill","trciHill")) {
      p[!ind] <- const + ptpareto(x[!ind],shape=1/EVTfit$gamma,scale=t,endpoint=EVTfit$endpoint) * (1-const)
      # CDF is 1 after endpoint
      p[x>EVTfit$endpoint] <- 1
    }
  
  }
  
  return(p)
}



# # Splicing quantiles
# SpliceQuant <- function(p, splicefit) {
#   
#   q <- numeric(length(p))
#   
#   MEfit <- splicefit$MEfit
#   EVTfit <- splicefit$EVTfit
#   
#   t <- splicefit$t
#   trunclower <- splicefit$trunclower
#   
#   const <- splicefit$const
#   
#   ind <- (p<const)
#   
#   # Quantiles of ME part
#   q[ind] <- ME_VaR(p[ind]/const, shape = MEfit$shape, alpha = MEfit$alpha, 
#                      theta = MEfit$theta, trunclower=trunclower, truncupper=t, 
#                    interval=c(trunclower,t)) 
#   # Quantiles of EVT part
#   if(splicefit$type=="GPD") {
#     q[!ind] <- t + EVTfit$sigma/EVTfit$gamma * ( ((1-p[!ind])/(1-const))^(-EVTfit$gamma) - 1 )
#   } else if(splicefit$type %in% c("Hill","cHill","ciHill")) {
#     q[!ind] <- t*((1-p[!ind])/(1-const))^(-EVTfit$gamma1)
#   } else if(splicefit$type %in% c("trHill","trciHill")) {
#     q[!ind] <- t*(1-(p[!ind]-const)/(1-const)*(1-(EVTfit$endpoint/t)^(-1/EVTfit$gamma)))^(-EVTfit$gamma)
#   }
#   
#   return(q)
# }


###########################################################################



# Plot of fitted survival function and ECDF estimator + bounds
SpliceECDF <- function(x, X, splicefit, alpha = 0.05, ...) {
  
  plot(x,1-SpliceCDF(x,splicefit=splicefit),type="l",xlab="x",ylab="1-F(x)", ...)

  # ECDF estimator
  fit  <- ecdf(X)
  # Empirical surivival function
  est <- 1-fit(x)
  
  # Confidence bounds using Dvoretzky-Kiefer-Wolfowitz inequality
  # http://stats.stackexchange.com/questions/55500/confidence-intervals-for-empirical-cdf
  n <- length(X)
  eps <- sqrt(1/(2*n)*log(2/alpha))
  lines(x,est,lty=1,col="red")
  lines(x,pmax(est-eps,0),lty=2,col="blue")
  lines(x,pmin(est+eps,1),lty=2,col="blue")
  legend("topright",c("Fitted survival function","Empirical survival function","95% confidence bounds"),
         lty=c(1,1,2),col=c("black","red","blue"))
}

# Plot of fitted survival function and Turnbull estimator + bounds
SpliceTB <- function(x, Z, I = Z, censored, splicefit, alpha = 0.05, ...) {
  
  plot(x,1-SpliceCDF(x,splicefit=splicefit),type="l",xlab="x",ylab="1-F(x)", ...)
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  
  
  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(Z==I,"right","interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain", conf.int=1-alpha)
  f = stepfun(fit$time,c(1,fit$surv))
  est <- f(x)
  
  lines(x,est,col="red")
  
  # Confidence bounds
  lines(fit$time,fit$lower, col = "blue", lty = 2)
  lines(fit$time,fit$upper, col = "blue", lty = 2)
  lines(x,est,lty=1,col="red")
  legend("topright",c("Fitted survival function","Turnbull estimator","95% confidence bounds"),
         lty=c(1,1,2),col=c("black","red","blue"))
}


# Probability - probability plot with ECDF
SplicePP <- function(x = sort(X), X, splicefit, log = FALSE, ...) {
  
  # ECDF estimator
  fit  <- ecdf(X)
  est <- 1-fit(x)
  
  if(log) {
    ind <- est>0
    plot(-log(est[ind]),-log(1-SpliceCDF(x[ind],splicefit=splicefit)),type="l",xlab="-log(Empirical survival probability)",
         ylab="-log(Fitted survival probability)", ...)
  } else {
    plot(est,1-SpliceCDF(x,splicefit=splicefit),type="l",xlab="Empirical survival probability",
         ylab="Fitted survival probability", ...)
  }
  abline(a=0,b=1)
}

# Probability - probability plot with Turnbull estimator
SplicePP_TB <- function(x = sort(Z), Z, I = Z, censored, splicefit, log=FALSE, ...) {
  
  # Sort the data with index return
  s <- sort(Z, index.return = TRUE)
  L <- s$x
  sortix <- s$ix
  R <- I[sortix]
  

  # Turnbull estimator
  event <- !censored[sortix]
  # 3 for interval censoring and 0 for right censoring
  event[event==0 & Z!=I] <- 3
  type <- ifelse(Z==I,"right","interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  f = stepfun(fit$time,c(1,fit$surv))
  est <- f(x)
  
  
  if(log) {
    plot(-log(est),-log(1-SpliceCDF(x,splicefit=splicefit)),type="l",xlab="-log(Empirical surv. probs)",
         ylab="-log(Fitted surv. probs)", ...)
  } else {
    plot(est,1-SpliceCDF(x,splicefit=splicefit),type="l",
         xlab="Empirical survival probability",ylab="Fitted probability", ...)
  }
  abline(a=0,b=1)
}


# Log-log plot with empirical survival function and fitted survival function
SpliceLL <- function(x = sort(X), X, splicefit, ...) {
  
  # ECDF estimator
  fit  <- ecdf(X)

  X <- sort(X)
  plot(log(X),log(1-fit(X)),ylab="log(emp. surv.)",xlab="log(X)",type="p", ...)
  lines(log(x),log(1-SpliceCDF(x,splicefit=splicefit)))
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
  type <- ifelse(Z==I,"right","interval")
  fit  <- survfit(Surv(time=L, time2=R, event=event, type=type) ~1, conf.type="plain")
  f = stepfun(fit$time,c(1,fit$surv))

  Z <- sort(Z)
  plot(log(Z),log(f(Z)),ylab="log(Turnbull surv.)",xlab="log(X)",type="p", ...)
  lines(log(x),log(1-SpliceCDF(x,splicefit=splicefit)))
}

