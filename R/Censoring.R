


# Hill estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary Hill estimator
cHill <- function(data, censored, plot = FALSE, add = FALSE, main = "Hill estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  censored <- checkCensored(censored, length(data))
  
  if(length(censored)!=1) {
    if(length(data) != length(censored)) {
      stop("data and censored should have the same length.")
    }
  } else {
    censored <- rep(0,length(data))
  }
  
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[s$ix])
  n <- length(X)
  Hill <- numeric(n)
  
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  K <- 1:(n-1)
  
  # Hill estimates
  # Fast vectorised version
  Hill[K] <- (cumsum(log(X[n-K+1]))- K*log(X[n-K])) / cumsum(delta[n-K+1]) 
  # Problems with division by 0
  Hill[K[cumsum(delta[n-K+1])==0]] <- NA
  

  # plots if TRUE
  plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  
  # output list with values of k and corresponding Hill estimates
  
  output(list(k=K, gamma1=Hill[K]),plot=plot,add=add)
  
}


# generalised Hill estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary generalised Hill estimator.
cgenHill <- function(data, censored, plot = FALSE, add = FALSE,  
         main = "Generalised Hill estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  censored <- checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  UH.scores <- numeric(n)
  Hill <- numeric(n)
  ind <- 1:(n-1)
  
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  gamma1 <- numeric(n)
  gamma1[ind] <- (cumsum(log(X[n-ind+1]))- ind*log(X[n-ind])) / cumsum(delta[n-ind+1]) 
  # Problems with division by 0
  indNA <- ind[cumsum(delta[n-ind+1])==0]
  gamma1[indNA] <- NA
  
  ### Generalised Hill estimates
  
  # Fast vectorised version
  UH.scores[ind] <- X[n-ind] * gamma1[ind]
  #Set to 1 to have no influence (we take logs!)
  UH.scores[indNA] <- 1
  
  
  # Stop at n-2 since UH.scores[n]=0, so log(UH.scores[n])=Inf
  K <- 1:max((n-2),1) 
  Hill[K] <- (cumsum(log(UH.scores[K]))- K*log(UH.scores[K+1])) / cumsum(delta[n-K+1]) 
  # Problems with division by 0
  Hill[K[cumsum(delta[n-K+1])==0]] <- NA
  

  # plots if TRUE
  plotfun(K, Hill[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  
  # output list with values of k and
  # corresponding Zipf estimates
  
  output(list(k=K, gamma1=Hill[K]), plot=plot, add=add)
  
}

# Moment estimator for censored data
# censored is a vector which is 1 if a data point is censored and 0 otherwise,
# giving censored=0 results in the ordinary Moment estimator
cMoment <- function(data, censored, plot = FALSE, add = FALSE, main = "Moment estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  censored <- censored <- checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[s$ix])
  n <- length(X)
  M1 <- numeric(n)
  M2 <- numeric(n)
  Mom <- numeric(n)
  K <- 1:(n-1)
  
  
  
  ### Moment estimates
  
  ######################
  #Fast vectorised version
  
  M1[K] <- cumsum(log(X[n-K+1]))/K - log(X[n-K])
  # for k=1:(n-1)
  # 1/k*sum_{j=1}^k [log(X[n-j+1])-log(X[n-k])]^2
  # = 1/k*sum_{j=1}^k log(X[n-j+1])^2 - 2*log(X[n-k])*1/k*sum_{j=1}^k log(X[n-j+1]) + log(X[n-k])^2
  M2[K] <- cumsum(log(X[n-K+1])^2)/K-2*cumsum(log(X[n-K+1]))/K*log(X[n-K]) + log(X[n-K])^2
  
  Mom <- M1 + 1 - (1-M1^2/M2)^(-1)/2
  #Add censoring
  Mom[K] <- Mom[K] / (cumsum(delta[n-K+1])/K )

  # Problems with division by 0
  Mom[K[cumsum(delta[n-K+1])==0]] <- NA
  
  ######################
  
  # plots if TRUE
  plotfun(K, Mom[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  
  # output list with values of k and
  # corresponding estimates for gamma1, b and beta
  output(list(k=K, gamma1=Mom[K]),plot=plot,add=add)
  
}


cGPDmle <- function(data, censored, start = c(0.1,1), warnings = FALSE, 
                 plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  censored <- checkCensored(censored, length(data))
  
  s <- sort(data, index.return = TRUE)
  X <- s$x
  sortix <- s$ix
  # delta is 1 if a data point is not censored, we also sort it with X
  delta <- !(censored[sortix]*1)
  n <- length(X)
  K <- 1:(n-1)
  
  # Proportion of non-censored
  pk <- cumsum(delta[n-K+1])/K
  
  # Apply ordinary POT estimator to the data
  L <- POT(data=X, start=start, warnings=warnings, plot=FALSE, add=FALSE)
  
  
  # Estimates
  gamma1 <- L$gamma[K]
  sigma1 <- L$sigma[K]
  
  # Censored version
  gamma1 <- gamma1 / pk
    
  # plots if TRUE
  plotfun(K, gamma1[K], type="l", xlab="k", ylab="gamma1", main=main, plot=plot, add=add, ...)
  
  output(list(k=K, gamma1=gamma1[K], sigma1=sigma1[K]), plot=plot, add=add)
}

cPOT <- cGPDmle

