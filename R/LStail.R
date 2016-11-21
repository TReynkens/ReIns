
# Computes the least squared estimates of gamma, b and rho 
# based on (5.36) in Beirlant et al. (2004) for a numeric vector of 
# observations (data) and as a function of k
#
# lambda chosen between 0 and 1

# Set rho or use estimator from Beirlant et al. (2002) when rho
# is set to NULL. We then need a value for lambda (between 0 and 1)

LStail <- function(data, rho = -1, lambda = 0.5, logk = FALSE, plot = FALSE, add = FALSE, 
                       main = "LS estimates of EVI", ...) {
  
  # Check input arguments
  .checkInput(data)
  
  X <- as.numeric(sort(data))
  n <- length(X)
  par <- matrix(nrow=n, ncol=3)
  UH.scores <- numeric(n)
  Hill <- numeric(n)
  K <- 1:(n-1)

  if (!is.numeric(rho) & !is.null(rho)) {
    stop("rho should be NULL or numeric.")
  }
  
  if (is.null(rho) & !is.numeric(lambda)) {
    stop("lambda should be numeric.")
  }
  
  
  
  if (!is.null(rho) & length(rho)!=1) {
    rho <- rho[1]
    warning("Input argument rho should be a numeric of length 1 or NULL, only first element of rho used.")
  }
  
  if (length(lambda)!=1) {
    lambda <- lambda[1]
    warning("Input argument lambda should be a numeric of length 1, only first element of lambda used.")
  }
  
  if (lambda<=0 | lambda >=1) {
    warning("Input argument lambda should be strictly between 0 and 1.")
  }
  
  ### LS estimates
  
  # Hill estimator for positive EVI
  # Fast vectorised version
  H <- numeric(n)
  H[K] <- cumsum(log(X[n-K+1])) / K - log(X[n-K])
  
  #UH scores (generalised QQ-plot)
  UH.scores[K] <- X[n-K] * H[K]
  
  
  # Estimate rho using single value or using estimator of Beirlant et al. (2002) 
  # when rho=NULL
  if (is.null(rho)) {
    
    # Generalised Hill estimator
    # Note that the generalised Hill estimator is LS estimator for model Z_j = gamma + eps_j
    # Stop at n-2 since UH.scores[n]=0, so log(UH.scores[n])=Inf
    ind <- 1:max((n-2), 1)
    Hill[ind] <- cumsum(log(UH.scores[ind])) / ind - log(UH.scores[ind+1])
    
    eps = sqrt(.Machine$double.eps)
    
    for (k in (n-1):1) {
      if (floor(lambda^2*k)>0 & floor(lambda^2*k)!=floor(lambda*k)) {
        h0 <- Hill[k]
        h1 <- Hill[floor(lambda*k)]
        h2 <- Hill[floor(lambda^2*k)]

        if (abs(h1-h2)>eps & abs(h1-h0)>eps & (h2-h1)*(h1-h0)>0 &
            is.finite(h1) & is.finite(h0)) {
          par[k, 3] <- -log((h2-h1)/(h1-h0))/log(lambda)
        }
      }
      
    }
  } else {
    par[,3] <- rho
  }

  # Z-values
  Z <- (K+1)*log(UH.scores[K]/UH.scores[K+1])
  
  ######################
  # (semi)-vectorised version
  
  # Obtain LS estimators for gamma and b using full model Z_j = gamma + b*(j/k)^(-rho) + eps_j
  if (is.null(rho)) {
    # Semi-vectorised
    for (k in (n-1):1) {
      #1/k*sum_{j=1}^k (j/k)^(-rho(k))* Z[j]
      par[k, 2] <- sum(((1:k)/k)^(-par[k, 3])*Z[1:k]) / k 
    }
    par[K, 2] <- ((1-par[K, 3])^2*(1-2*par[K, 3])/(par[K, 3])^2)*(par[K, 2]-1/(1-par[K, 3])*cumsum(Z)/K)
  } else {
    # Vectorised since par[K,3] is constant (rho)
    # for k=1:(n-1)
    # 1/k*sum_{j=1}^k [(j/k)^(-rho(k))-1/(1-rho)] * Z[j]
    # = k^(rho-1)*sum_{j=1}^k j^(-rho)*Z[j] - 1/(1-rho)*1/k*sum_{j=1}^k Z[j] 
    par[K, 2] <- (1-rho)^2 * (1-2*rho) / rho^2 * (cumsum(K^(-rho)*Z) * K^(rho-1) - 1/(1-rho) * cumsum(Z) / K)
  }
    
  par[K, 1] <- cumsum(Z) / K - par[K, 2] / (1-par[K, 3])
  

  
  ######################
  if (logk) {
    .plotfun(log(K), par[K, 1], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, par[K, 1], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }
  
  # output list with values of k and
  # corresponding estimates for gamma, b and rho
  
  .output(list(k=K, gamma=par[K, 1], b=par[K, 2], rho=par[K, 3]), plot=plot, add=add)

}


TSfraction <- LStail