
###############################################################
# Pareto

dpareto <- function(x, shape, scale = 1, log = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  d <- ifelse(x>=scale, shape/scale*(scale/x)^(shape+1), 0)
  
  if (log) d <- log(d)
  
  return(d)
  
}

ppareto <- function(x, shape, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  p <- ifelse(x>=scale, 1-(scale/x)^shape, 0)
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  
  return(p)
  
}

qpareto <- function(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (all(p>=0 & p<=1)) {
    
    return(scale*(1-p)^(-1/shape))
    
  } else {
    stop("p should be between 0 and 1.")
  }
  
}

rpareto <- function(n, shape, scale = 1) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  return(qpareto(runif(n),scale=scale,shape=shape))
}


###############################################################
# Truncated Pareto

dtpareto <- function(x, shape, scale=1, endpoint=Inf, log = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  if(endpoint<=scale) {
    stop("endpoint should be strictly larger than scale.")
  }
  
  d <- ifelse(x<=endpoint, 
              dpareto(x,shape=shape,scale=scale)/ppareto(endpoint,shape=shape,scale=scale),
              0)
  if (log) d <- log(d)
  
  return(d)
}

ptpareto <- function(x, shape, scale=1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  if(endpoint<=scale) {
    stop("endpoint should be strictly larger than scale.")
  }
 
  p <- ifelse(x<=endpoint, 
              ppareto(x,shape=shape,scale=scale)/ppareto(endpoint,shape=shape,scale=scale),
              1)
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
  
  #return(ifelse(x>=scale,(1-(scale/x)^shape)/(1-(scale/endpoint)^shape),0))
}

qtpareto <- function(p, shape, scale=1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  if(endpoint<=scale) {
    stop("endpoint should be strictly larger than scale.")
  }
  
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (all(p>=0 & p<=1)) {

    return(qpareto(p*ppareto(endpoint,shape=shape,scale=scale),shape=shape,scale=scale))
    #return(scale*(1-p*(1-(scale/endpoint)^shape))^(-1/shape))
  } else {
    stop("p should be between 0 and 1.")
  }
  
}

rtpareto <- function(n, shape, scale=1, endpoint=Inf) {
  if (shape<=0) {
    stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
    stop("scale should be strictly positive.")
  }
  
  if(endpoint<=scale) {
    stop("endpoint should be strictly larger than scale.")
  }
  
  return(qtpareto(runif(n),scale=scale,shape=shape,endpoint=endpoint))
}


######################
# Generalised Pareto distribution (GPD)

eps <- 10^(-14)


dgpd <- function(x, gamma, mu = 0, sigma, log = FALSE) {
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  if (abs(gamma)<eps) {
    d <- ifelse(x>=mu, 1/sigma * exp(-(x-mu)/sigma), 0)
  } else {
    d <- ifelse(x>=mu, 1/sigma * (1+gamma*(x-mu)/sigma)^(-1/gamma-1), 0)
  }
 
  if(gamma < -eps) {
    d[x>mu-sigma/gamma] <- 0
  } 
  
  if (log) d <- log(d)
  
  return(d)
}

pgpd <- function(x, gamma, mu = 0, sigma, lower.tail = TRUE, log.p = FALSE) {
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  if (abs(gamma)<eps) {
    p <- ifelse(x>=mu, 1-exp(-(x-mu)/sigma), 0)
  } else {
    p <- ifelse(x>=mu, 1-(1+gamma*(x-mu)/sigma)^(-1/gamma), 0)
  }
  
  if(gamma < -eps) {
    p[x>=mu-sigma/gamma] <- 1
  }
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  
  return(p)
}


qgpd <- function(p, gamma, mu = 0, sigma, lower.tail = TRUE, log.p = FALSE) {
  
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (!all(p>=0 & p<=1)) {
    stop("p should be between 0 and 1.")
  }
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  

  if (abs(gamma)<eps) {
    q <- mu - sigma * log(1-p)
  } else {
    q <- mu + sigma/gamma * ((1-p)^(-gamma) - 1)
  }
  
  return(q)
}


rgpd <- function(n, gamma, mu = 0, sigma) {
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  qgpd(runif(n), gamma=gamma, mu=mu, sigma=sigma)
}

######################
# Truncated GPD

dtgpd <- function(x, gamma, mu = 0, sigma, endpoint=Inf, log = FALSE) {
  
  d <- ifelse(x<=endpoint, 
              dgpd(x,gamma=gamma,mu=mu,sigma=sigma)/pgpd(endpoint,gamma=gamma,mu=mu,sigma=sigma),
              0)
  
  if (log) d <- log(d)
  
  return(d)
  
}

ptgpd <- function(x, gamma, mu = 0, sigma, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  if(endpoint<=mu) {
    stop("endpoint should be strictly larger than mu.")
  }
  
  p <- ifelse(x<=endpoint, 
              pgpd(x,gamma=gamma,mu=mu,sigma=sigma)/pgpd(endpoint,gamma=gamma,mu=mu,sigma=sigma),
              1)
  
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  
  return(p)
}


qtgpd <- function(p, gamma, mu = 0, sigma, endpoint = Inf, lower.tail = TRUE, log.p = FALSE) {
  
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (!all(p>=0 & p<=1)) {
    stop("p should be between 0 and 1.")
  }
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  if(endpoint<=mu) {
    stop("endpoint should be strictly larger than mu.")
  }

  q <- qgpd(p*pgpd(endpoint,gamma=gamma,mu=mu,sigma=sigma),gamma=gamma,mu=mu,sigma=sigma)
  
  return(q)
   
}


rtgpd <- function(n, gamma, mu = 0, sigma, endpoint = Inf) {
  
  if (sigma<0) {
    stop("sigma should be strictly positive.")
  }
  
  if(endpoint<=mu) {
    stop("endpoint should be strictly larger than mu.")
  }
  
  
  qtgpd(runif(n),gamma=gamma,mu=mu,sigma=sigma, endpoint=endpoint)
}


###############################################################
# Burr (type XII)

dburr = function(x, alpha, rho, log = FALSE) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  d <- ifelse(x>0, alpha * (1+x^(-rho*alpha))^(1/rho-1) * x^(-rho*alpha-1), 0)
  
  if (log) d <- log(d)
  
  return(d)
  
}


pburr = function(x, alpha, rho, lower.tail = TRUE, log.p = FALSE) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  p <- ifelse(x>0, 1-(1+x^(-rho*alpha))^(1/rho), 0)
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
  
}

qburr = function(p, alpha, rho, lower.tail = TRUE, log.p = FALSE) {
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  if (all(p>=0 & p<=1)) {
    return( ((1-p)^rho-1)^(-1/(rho*alpha)) )
    
  } else {
    stop("p should be between 0 and 1.")
  }
  
}

rburr = function(n, alpha, rho) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  return(qburr(runif(n), rho=rho, alpha=alpha))
}


###############################################################
# Truncated Burr


dtburr = function(x, alpha, rho, endpoint=Inf, log = FALSE) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  d <- ifelse(x<=endpoint, dburr(x, alpha=alpha, rho=rho)/pburr(endpoint, alpha=alpha, rho=rho), 0)
  
  if (log) d <- log(d)
  
  return(d)
  
}


ptburr = function(x, alpha, rho, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  #   A = 1-(1+endpoint^(-rho*alpha))^(1/rho)
  #   return(ifelse(x>0 & x<endpoint,(1-(1+x^(-rho*alpha))^(1/rho))/A,0))
  
  p <- ifelse(x<=endpoint, pburr(x, alpha=alpha, rho=rho)/pburr(endpoint, alpha=alpha, rho=rho), 1)
  
  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
  
}


qtburr = function(p, alpha, rho, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  if (all(p>=0 & p<=1)) {
    #     A = 1-(1+endpoint^(-rho*alpha))^(1/rho)
    #     return(((1-p*A)^rho-1)^(-1/(rho*alpha)))
    return( qburr(p*pburr(endpoint, alpha=alpha, rho=rho), alpha=alpha, rho=rho) )
    
  } else {
    stop("p should be between 0 and 1.")
  }
  
}

rtburr = function(n, alpha, rho, endpoint=Inf) {
  
  if (alpha<=0) {
    stop("alpha should be strictly positive.")
  }
  
  if (rho>=0) {
    stop("rho should be strictly negative")
  }
  
  if(endpoint<=0) {
    stop("endpoint should be strictly larger than 0.")
  }
  
  return(qtburr(runif(n),alpha=alpha,rho=rho,endpoint=endpoint))
}

###############################################################
# Truncated log-normal

dtlnorm <- function(x, meanlog = 0, sdlog = 1, endpoint=Inf, log = FALSE) {
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  d <- ifelse(x<=endpoint, 
              dlnorm(x,meanlog=meanlog,sdlog=sdlog)/plnorm(endpoint,meanlog=meanlog,sdlog=sdlog), 
              0)
  
  if (log) d <- log(d)
  
  return(d)
  #return(dnorm((log(x)-meanlog)/sdlog)/pnorm((log(endpoint)-meanlog)/sdlog))
}

ptlnorm <- function(x, meanlog = 0, sdlog = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {

  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  p <- ifelse(x<=endpoint, 
              plnorm(x,meanlog=meanlog,sdlog=sdlog)/plnorm(endpoint,meanlog=meanlog,sdlog=sdlog), 
              1)

  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
  #ifelse(x<=0, 0, pnorm((log(x)-meanlog)/sdlog)/pnorm((log(endpoint)-meanlog)/sdlog))
}

qtlnorm <- function(p, meanlog = 0, sdlog = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE)  {
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  if (all(p>=0 & p<=1)) {

    return(qlnorm(p*plnorm(endpoint,meanlog=meanlog,sdlog=sdlog),meanlog=meanlog,sdlog=sdlog))
#     A = pnorm((log(endpoint)-meanlog)/sdlog)
#     return(exp(meanlog+sdlog*qnorm(p*A)))
  } else {
    stop("p should be between 0 and 1.")
  }
  
}

rtlnorm <- function(n, meanlog = 0, sdlog = 1, endpoint=Inf) {
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  return(qtlnorm(runif(n),meanlog=meanlog,sdlog=sdlog,endpoint=endpoint))
}

# ###############################################################
# # Weibull 
# 
# dweibull <- function(x, lambda, tau = 1, log = FALSE) {
#   
#   if (lambda<=0) {
#    stop("lambda should be strictly positive.")
#   }
#   
#   if (tau<=0) {
#    stop("tau should be strictly positive.")
#   }
#   
# 
#   d <- ifelse(x<=0, 0, lambda*tau*x^(tau-1)*exp(-lambda*x^tau))
#   
#   if (log) d <- log(d)
#   
#   return(d)
# }
# 
# 
# pweibull <- function(x, lambda, tau = 1, lower.tail = TRUE, log.p = FALSE) {
#   
#   if (lambda<=0) {
#    stop("lambda should be strictly positive.")
#   }
#   
#   if (tau<=0) {
#    stop("tau should be strictly positive.")
#   }
#   
#   
#   p <- ifelse(x<=0, 0, 1-exp(-lambda*x^tau))
#   
#   if (!lower.tail) p <- 1-p
#   
#   if (log.p) p <- log(p)
#   
#   return(p)
# }
# 
# 
# qweibull <- function(p, lambda, tau = 1, lower.tail = TRUE, log.p = FALSE) {
#   
#   if (lambda<=0) {
#    stop("lambda should be strictly positive.")
#   }
#   
#   if (tau<=0) {
#    stop("tau should be strictly positive.")
#   }
#   
#   if (log.p) p <- exp(p)
#   
#   if (!lower.tail) p <- 1-p
#   
#   
#   if (all(p>=0 & p<=1)) {
# 
#     return((-log(1-p)/lambda)^(1/tau))
#   } else {
#     stop("p should be between 0 and 1.")
#   }
# 
# }
# 
# rweibull <- function(n, lambda, tau = 1) {
#   
#   if (lambda<=0) {
#    stop("lambda should be strictly positive.")
#   }
#   
#   if (tau<=0) {
#    stop("tau should be strictly positive.")
#   }
#   
#   return(qweibull(runif(n),lambda=lambda,tau=tau))
#   
# }

###############################################################
#Truncated Weibull

dtweibull <- function(x, shape, scale = 1, endpoint=Inf, log = FALSE) {
  
  if (shape<=0) {
   stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
   stop("scale should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  d <- ifelse(x<=endpoint,
              dweibull(x,shape=shape,scale=scale)/pweibull(endpoint,shape=shape,scale=scale), 
              0)
  
  if (log) d <- log(d)
  
  return(d)
}

ptweibull <- function(x, shape, scale = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (shape<=0) {
   stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
   stop("scale should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  p <- ifelse(x<=endpoint, 
              pweibull(x,shape=shape,scale=scale)/pweibull(endpoint,shape=shape,scale=scale), 
              1)
  #ifelse(x<=0, 0, (1-exp(-(x/beta)^alpha))/1-exp(-(endpoint/beta)^alpha))

  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}

qtweibull <- function(p, shape, scale = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (shape<=0) {
   stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
   stop("scale should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  
  
  if (all(p>=0 & p<=1)) {

    return(qweibull(p*pweibull(endpoint,shape=shape,scale=scale),shape=shape,scale=scale))
  } else {
    stop("p should be between 0 and 1.")
  }
}

rtweibull <- function(n, shape, scale = 1, endpoint=Inf) {
  
  if (shape<=0) {
   stop("shape should be strictly positive.")
  }
  
  if (scale<=0) {
   stop("scale should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  return(qtweibull(runif(n),shape=shape,scale=scale,endpoint=endpoint))
}

###############################################################
#Truncated Exponential

dtexp <- function(x, rate = 1, endpoint=Inf, log = FALSE) {
  
  if (rate<=0) {
    stop("rate should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  d <- ifelse(x<=endpoint, dexp(x,rate=rate)/pexp(endpoint,rate=rate), 0)
  
  if (log) d <- log(d)
  
  return(d)
}


ptexp <- function(x, rate = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (rate<=0) {
    stop("rate should be strictly positive.")
  }
  
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  p <- ifelse(x<=endpoint, pexp(x,rate=rate)/pexp(endpoint,rate=rate), 1)

  if (!lower.tail) p <- 1-p
  
  if (log.p) p <- log(p)
  
  return(p)
}

qtexp <- function(p, rate = 1, endpoint=Inf, lower.tail = TRUE, log.p = FALSE) {
  
  if (rate<=0) {
    stop("rate should be strictly positive.")
  }
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  
  
  if (log.p) p <- exp(p)
  
  if (!lower.tail) p <- 1-p
  

  if (all(p>=0 & p<=1)) {
    
    return( qexp(p*pexp(endpoint,rate=rate),rate=rate))
    
  } else {
    stop("p should be between 0 and 1.")
  }
}

rtexp <- function(n, rate = 1, endpoint=Inf) {
  
  if (rate<=0) {
    stop("rate should be strictly positive.")
  }
  
  
  
  if(endpoint<=0) {
    stop("endpoint should be strictly positive.")
  }
  

  return(qtexp(runif(n),rate=rate,endpoint=endpoint))
}

