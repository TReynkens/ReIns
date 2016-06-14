
# This file contains implementations of estimators that are suitable for a truncated
# Pareto distribution.

###################################################################################################

# Estimation of extreme value index of truncated Pareto model
trHill <- function(data, r = 1, tol = 1e-8, maxiter = 100, logk = FALSE, plot = FALSE, add = FALSE, 
                         main="Estimates of EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  gamma <- numeric(n)
  H <- numeric(n)
  H1 <- numeric(n)
  R <- numeric(n)
  K <- r:(n-1)
  
  # Truncated Hill estimator 
  
  # Trimmed Hill
  H[K] <- cumsum(log(X[(n-r+1):(n-tail(K,1)+1)]))/(K-r+1) - log(X[n-K]) 
  K1 <- 1:(n-1)
  # Hill
  H1[K] <- (cumsum(log(X[n-K1+1]))/K1 - log(X[n-K1]))[K] 
  
  R[K] <- X[n-K]/X[n-r+1]
  
  # Newton-Raphson to obtain estimates of gamma
  for (k in (n-1):r) {  
    # Note that y = gamma = 1/alpha
    #y <- H1[k]
    #Use Hill's trimmed estimator as initial approximation
    y <- H[k]
    i <- 1
    while (i <= maxiter){
      #temp[i] <- y
      ym <- 1/y
      z <- y + (H[k]-y-R[k]^ym*log(R[k])/(1-R[k]^ym))/(1-ym^2*R[k]^ym*log(R[k])^2/((1-R[k]^ym)^2))
      
      if (abs(z-y) < tol | identical(abs(z-y),NaN)) { break }
      
      i <- i+1
      y <- z
    }
    gamma[k] <- z
  }
  
  gamma[gamma<=0] <- NA
  
  # plots if TRUE  
  if (logk) {
    .plotfun(log(K), gamma[K], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, gamma[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }
  
  # output list with values of k and
  # corresponding estimates for gamma
  
  .output(list(k=K, gamma=gamma[K], H=H[K]),plot=plot,add=add)
  
}


# Estimator for odds ratio DT
trDT <- function(data, r = 1, gamma, plot=FALSE, add=FALSE, 
                      main="Estimates of DT", ...) {
  
  # Check input arguments
  .checkInput(data,gamma=gamma,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  DT <- numeric(n)
  K <- r:(n-1)
  
  if(length(gamma)!=length(K)) {
    stop(paste("gamma should have length", length(K)))
  }
  
  
  A <- 1/gamma
  R <- X[n-K]/X[n-r+1]
  
  DT[K] <-  pmax( (K+1)/(n+1) * (R^A-r/(K+1)) / (1-R^A), 0)
  
  # plots if TRUE  
  .plotfun(K, DT[K], type="l", xlab="k", ylab=expression(D[T]), main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k and
  # corresponding estimates for DT
  
  .output(list(k=K, DT=DT[K]),plot=plot,add=add)
  
}

# Estimator for endpoint
trEndpoint <- function(data, r = 1, gamma, plot = FALSE, add = FALSE, 
                      main = "Estimates of Endpoint", ...) {
  
  # Check input arguments
  .checkInput(data,gamma=gamma,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  Tk <- numeric(n)
  K <- r:(n-1)
  

  
  if(length(gamma)!=length(K)) {
    stop(paste("gamma should have length", length(K)))
  }
  
  R <- X[n-K] / X[n]
  
  Tk[K] <- pmax(X[n-K] * ((R^(1/gamma)-1/(K+1)) / (1-1/(K+1)))^(-gamma), X[n])
  
  # plots if TRUE  
  .plotfun(K, Tk[K], type="l", xlab="k", ylab=expression(T[k]), main=main, plot=plot, add=add, ...)
  
  # output list with values of k and
  # corresponding estimates for DT
  
  .output(list(k=K, Tk=Tk[K]),plot=plot,add=add)
  
}

# Estimator for extremes quantile
# rough indicates if we are in the rough truncation case
trQuant <- function(data, r = 1, rough = TRUE, gamma, DT, p, plot = FALSE, add = FALSE,
                         main="Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data,gamma=gamma,DT=DT,r=r)
  .checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  quant <- numeric(n)
  K <- r:(n-1)
  
  # Estimator for extreme quantiles

  if (rough) {
    #quant[K] <- X[n-K] * ((K+1)/((n+1)*p))^gamma * ((1+DT*(n+1)/(K+1))/(1+DT/p))^gamma
    quant[K] <- X[n-K] *  ((DT+(K+1)/(n+1))/(DT+p))^gamma
    
  } else {
    quant[K] <- X[n-K] * ((K+1)/(p*(n+1)))^gamma
  }


  
  # plots if TRUE
  
  .plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p),plot=plot,add=add)
  
}

# Estimator for extreme quantile of original data W
trQuantW <- function(data, gamma, DT, p, plot = FALSE, add = FALSE,
                    main="Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data,gamma=gamma,DT=DT)
  .checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  quant <- numeric(n)
  K <- 1:(n-1)
  
  # Estimator for extreme quantiles
  
  
  quant[K] <- X[n-K] * ((DT+(K+1)/(n+1))/(p*(1+DT)))^gamma[K]
    
 
  
  # plots if TRUE
  
  .plotfun(K, quant[K], type="l", xlab="k", ylab=expression(Q[W](1-p)), main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=quant[K], p=p),plot=plot,add=add)
  
}




# Estimator of exceedance probability
trProb <- function(data, r = 1, gamma, q, warnings = TRUE, plot = FALSE, add = FALSE,
                         main="Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data,gamma=gamma,r=r)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  prob <- numeric(n)
  K <- r:(n-1)
  
  R <- X[n-K]/X[n-r+1]
  
  # Estimator for exceedance probability
  prob[K] <- (K+1)/(n+1) * ( (q/X[n-K])^(-1/gamma) - R^(1/gamma) ) / (1-R^(1/gamma))
  
  if (q > max(data) & warnings) {
    warning("q should be smaller than the largest observation.")
  } 

  
  prob[prob<0 | prob>1] <- NA
  
  # plots if TRUE
  
  if( !all(is.na(prob[K])) ) {
    .plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  }
  
  
  # output list with values of k, exceedance probabilitye estimates 
  # and the considered high quantile q
  
  .output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
  
}

# Pareto QQ-plot adapted for truncation
trParetoQQ <- function(data, r = 1, DT, kstar = NULL, plot = TRUE, main = "TPa QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(data, DT=DT)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  j <- 1:n
  K <- r:(n-1) 
  
  if(length(DT)==1) {
    DT <- rep(DT,length(K))
  }
  
  # calculate theoretical and empirical quantiles
  pqq.emp <- log(X[n-j+1])
  
  # Select kstar based on maximising correlation
  if(is.null(kstar)) {
    
    # Ignore first values for k if K is large
    k <- ifelse(length(K)>10, 11, 1):length(K)

    # Compute correlation between log(X[n-j+1]) and DT[i]+j/n
    cors <- numeric(length(k)) 
    
    for(i in 1:length(k)) {
      cors[i] <- abs( cor(log(X[n-(1:k[i])+1]), -log(DT[K==k[i]]+(1:k[i])/n)) )
    }
    # Value for k which maximises correlation
    kstar <- k[which.max(cors)]
    
  } else {
    
    # Check if input for kstar is valid
    if (!is.numeric(kstar) | length(kstar)>1) {
      stop("kstar should be a numeric of length 1.")
    }
    
    if (kstar<=0 | !.is.wholenumber(kstar)) {
      stop("kstar should be a strictly positive integer.")
    }
    
    if (kstar>n-1) {
      stop(paste0("kstar should be strictly smaller than ",n,"."))
    }
    
    if (kstar<=10) {
      warning("kstar should be strictly larger than 10.")
    }
    
  }
  
  pqq.the <- -log(DT[K==kstar]+j/(n+1))

  # Old margins
  oldmar <-  par("mar")
  
  # Increase margins to the left
  par(mar=c(5.1, 4.6, 4.1, 2.1))
  
  # plots if TRUE
  xlab <- bquote(-log(hat(D)[.(paste0("T,",r,",",kstar))]+j/(n+1)))
  ylab <- bquote(log(X["n-j+1,n"]))
  .plotfun(pqq.the, pqq.emp, type="p", xlab=xlab, ylab=ylab, 
          main=main, plot=plot, add=FALSE, ...)
  
  # Reset margins
  par(mar=oldmar)
  
  
  .output(list(pqq.the=pqq.the, pqq.emp=pqq.emp, kstar=kstar, DTstar = DT[K==kstar]), plot=plot, add=FALSE)
}


# Test for non-truncated vs. truncated Pareto-type tails
trTest <- function(data, alpha = 0.05, plot = TRUE, main = "Test for truncation", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  
  K <- 1:(n-1)
  
  H <- Hill(X)$gamma
  
  # Compute L(H_{k,n})
  E <- numeric(n)

  for(k in K) {
    E[k] <- sum( (X[n-k]/X[n-(1:k)+1])^(1/H[k]) ) / k
  }  
    
  L <- (E[K]-0.5)/(1-E[K])
  
  # Test value
  tv <- sqrt(12*K)*L
  
  # Critical value
  cv <- qnorm(alpha)
  
  # Reject
  reject <- (tv<=cv)
  
  # P-values
  Pval <- pnorm(tv)
  
  # Plot P-values 
  .plotfun(K, Pval[K], ylim=c(0,1), type="l", xlab="k", ylab="P-value", main=main, plot=plot, add=FALSE, ...)
  if (plot) abline(h=alpha, col="blue", ...)
  
  .output(list(k=K, testVal=tv, critVal=cv, Pval=Pval, reject=reject), 
         plot=plot, add=FALSE)
  
}
