
# Computes the Hill estimates of gamma (Section 4.2) 
# for a numeric vector of observations (data) and as a
# function of k
#
# K indicates if the estimates are plotted as a function of k or log(X[n-k])
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot

Hill <- function(data, k = TRUE, plot = FALSE, add = FALSE, main = "Hill estimates of EVI", ...) {
	
  # Check input arguments
  .checkInput(data)
  
  n <- length(data)
  Hill <- numeric(n)
  X <- as.numeric(sort(data))
  
	K <- 1:(n-1)
	
  # Hill estimates
  
	######################
  # Vectorised version
  Hill[K] <- cumsum(log(X[n-K+1])) / K - log(X[n-K])
  
  ######################
  # plots if TRUE
  
  if (k) {
    .plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(log(X[n-K]), Hill[K], type="l", xlab=bquote(log(X["n-k,n"])), 
            ylab="gamma", main=main, plot=plot, add=add, ...)
  }
  
  
  # output list with values of k and corresponding Hill estimates
  
  .output(list(k=K, gamma=Hill[K]),plot=plot,add=add)

}

##########################################################################

# Computes estimates of small exceedance probability 1-F(q) 
# (Section 4.6.1) for a numeric vector of observations 
# (data) and as a function of k
#
# Estimates are based on prior estimates gamma for the
# extreme value index (Hill)
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot

Prob <- function(data, gamma, q, plot = FALSE, add = FALSE, 
                 main = "Estimates of small exceedance probability", ...) {
  
  # Check input arguments
  .checkInput(data,gamma)
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  wp <- numeric(n)
  K <- 1:(n-1)
  
  # Weissman estimator for probabilities
  
  wp[K] <- (K+1)/(n+1) * (q/X[n-K])^(-1/gamma[K])
  wp[wp<0 | wp>1] <- NA
  
  # plots if TRUE
  .plotfun(K, wp[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, P=wp[K], q=q),plot=plot,add=add)
}

Weissman.p <- Prob


# Return period 
Return <- function(data, gamma, q, plot = FALSE, add = FALSE, 
                   main = "Estimates of return period", ...) {
  # Check input arguments
  .checkInput(data,gamma)
  
  if (length(q)>1) {
    stop("q should be a numeric of length 1.")
  }
  
  X <- as.numeric(sort(data))
  n <- length(X)
  wr <- numeric(n)
  K <- 1:(n-1)
  
  # Weissman estimator for return period 
  
  wr[K] <- (n+1)/(K+1) * (q/X[n-K])^(1/gamma[K])
  wr[wr<1] <- NA
  
  # plots if TRUE
  .plotfun(K, wr[K], type="l", xlab="k", ylab="1/(1-F(x))", main=main, plot=plot, add=add, ...)
  
  
  # output list with values of k, corresponding return period estimates 
  # and the considered large quantile q
  
  .output(list(k=K, R=wr[K], q=q),plot=plot,add=add)
}

Weissman.r <- Return 

##########################################################################

# Computes estimates of extreme quantile Q(1-p) 
# (Section 4.6.1) for a numeric vector of observations 
# (data) and as a function of k
#
# Estimates are based on prior estimates gamma for the
# extreme value index (Hill)
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot

Quant <- function(data, gamma, p, plot = FALSE, add = FALSE, 
                  main = "Estimates of extreme quantile", ...) {
  
  # Check input arguments
  .checkInput(data,gamma)
  
  .checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  wq <- numeric(n)
  K <- 1:(n-1)
  
  # Weisman estimator for quantiles
  
  wq[K] <- X[n-K] * ((K+1)/((n+1)*p))^(gamma[K])
  
  # plots if TRUE
  .plotfun(K, wq[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  # output list with values of k, corresponding quantile estimates 
  # and the considered small tail probability p
  
  .output(list(k=K, Q=wq[K], p=p),plot=plot,add=add)
}

Weissman.q <- Quant