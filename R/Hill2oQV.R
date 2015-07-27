# Computes bias reduced estimates of gamma based on the
# quantile view (Section 4.5.1) for a numeric vector of 
# observations (data) and as a function of k
#
# Method is based on an exponential regression model for
# log-spacings (Section 4.4) and maximum likelihood estimation
# on the parameters gamma, b and beta
#
# start contains the starting values for the ml optimization
# routine and can be altered if necessary
#
# Subsequent estimates are constrained by certain smoothness
# constraints on the parameter space
#
# If plot=TRUE then the estimates of gamma are plotted as a
# function of k
#
# If add=TRUE then the estimates of gamma are added to an existing
# plot


Hill.2oQV <- function(data, start = c(1,1,1), warnings = FALSE, plot = FALSE, add = FALSE, 
                      main = "Estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  H <- matrix(nrow=n, ncol=3) #empty matrix
  H[n,] <- start
  K <- 1:(n-1)

  
  # Hill 2nd order, quantile view
  mif <- 1.1					
  rhoh.min <- 0.5			
  
  
  #Numerical constant
  eps <- .Machine$double.eps
  
  
  # calculate minus log-likelihood of exponential and its derivative
  
  llsum <- function(w, y, x) {	
    # Zj ~  gamma + b_n,k * (j/(k+1))^beta
    lambda <- w[1] + w[2]*x^w[3]	
    if (min(lambda)< eps) lambda <- 10^10 #high value!
    sum(y/lambda + log(lambda))
  }
  
  llsum.dif <- function(w, y, x) {
    lambda <- w[1] + w[2]*x^w[3]
    
    if (min(lambda)<0) {
      part1 <- part2 <- part3 <- 0
    } else {
      # partial derivatives
      part1 <- sum(1/lambda - y/lambda^2) # wrt to gamma
      part2 <- sum(x^w[3]/lambda - (y*x^w[3])/lambda^2) # wrt b
      part3 <- sum((w[2]*log(x)*x^w[3])/lambda - (y*w[2]*log(x)*x^w[3])/lambda^2) # wrt to beta
    }
    
    return(c(part1, part2, part3))
  }
  
  Z <- (1:(n-1))*(log(X[n:2])-log(X[(n-1):1])) 
  
  
  for (k in (n-1):1) {
    xval <- (1:k)/(k+1)
    Zval <- Z[1:k]
    
    if (X[n-k]>0) { 
      # minimise llsum, give him the gradient, method, data ... 
      estim <- optim(par=start, fn=llsum, gr=llsum.dif, method = "L-BFGS-B", y=Zval, x=xval, 
                     lower=c(0.001, -mif*abs(H[k+1,2]), rhoh.min), 
                     upper=c(Inf, mif*abs(H[k+1,2]), mif*H[k+1,3]))
      if (estim$convergence>0 & warnings) {
        warning("Optimisation did not complete succesfully.")
        if (!is.null(estim$message)) {
          print(estim$message)
        }
      }
      estimators <- as.vector(estim$par)
      # Use estimated parameters as starting values for next value of k
      start <- estimators 
      
    } else {
      estimators <- rep(NA,3)
    } 
    H[k,] <- estimators
  }
  
  # plots if TRUE
  plotfun(K, H[K,1], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, gamma=H[K,1], b=H[K,2], beta=H[K,3]),plot=plot,add=add)

}




##########################################################################

# Computes second order refined estimates of extreme 
# quantile Q(1-p) (Section 4.6.2) based on the
# quantile view for a numeric vector of 
# observations (data) and as a function of k
#
# Method is based on a prior exponential regression model for
# log-spacings (Section 4.4) and maximum likelihood estimation
# on the parameters gamma, b and beta (Hill.2oQV)
#
# If plot=TRUE then the estimates of gamma are plotted as a
# function of k
#
# If add=TRUE then the estimates of gamma are added to an existing
# plot

Quant.2oQV <- function(data, gamma, b, beta, p, plot = FALSE, add = FALSE, 
                       main="Estimates of extreme quantile", ...) {
  
  # Check input arguments
  checkInput(data,gamma)
  
  checkProb(p)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  wq <- numeric(n)
  K <- 1:(n-1)
  
  # Weissman 2nd order estimator for quantiles, quantile view
  
  wq[K] <- X[n-K] * ((K+1)/((n+1)*p))^(gamma[K]) * exp( b[K]*(1-((K+1)/((n+1)*p))^(-beta[K]))/beta[K] )
  
  ### plots if TRUE
  
  # plots if TRUE
  plotfun(K, wq[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
  
  ### output list with values of k, corresponding quantile estimates 
  ### and the considered small tail probability p
  output(list(k=K, Q=wq[K], p=p),plot=plot,add=add)
}

Weissman.q.2oQV <- Quant.2oQV 

