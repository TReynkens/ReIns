

############################################################################################################

# Computes the moment estimates of gamma (Section 5.2.2) 
# for a numeric vector of observations (data) 
# and as a function of k
#
# If plot=TRUE then the estimates are plotted as a
# function of k
#
# If add=TRUE then the estimates are added to an existing
# plot


Moment <- function(data, plot = FALSE, add = FALSE, main = "Moment estimates of EVI", ...) {
  
  # Check input arguments
  checkInput(data)
  
  X <- as.numeric(sort(data))
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
  

  
  ######################
  
  # plots if TRUE
  plotfun(K, Mom[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  
  # output list with values of k and
  # corresponding estimates for gamma, b and beta
  output(list(k=K, gamma=Mom[K]),plot=plot,add=add)
  
  
}

#Same expressions as using the generalised Hill estimator
QuantMOM <- QuantGH
ProbMOM <- ProbGH
ReturnMOM <- ReturnGH


# #Estimate extreme quantiles using the moment estimator for gamma
# QuantMOM <- function(data, gamma, p, plot = FALSE, add = FALSE, 
#                       main = "Estimates of extreme quantile", ...) {
#   
#   # Check input arguments
#   checkInput(data,gamma,gammapos=FALSE)
#   checkProb(p)
#   
#   
#   X <- sort(data)
#   n <- length(X)
#   a <- numeric(n)
#   quant <- numeric(n)
#   K <- 1:(n-1)
# 
#   # Estimator for extreme quantiles
#   
#   H <- Hill(X)$gamma
#   a <- X[n-K]*H[K]*(1-pmin(gamma,0))
#   
#   quant[K] <- X[n-K] + a/gamma[K] * ( ((K+1)/((n+1)*p))^gamma[K] - 1 )
#     
#   # plots if TRUE
#   plotfun(K, quant[K], type="l", xlab="k", ylab="Q(1-p)", main=main, plot=plot, add=add, ...)
#   
#   # output list with values of k, corresponding quantile estimates 
#   # and the considered small tail probability p
#   
#   output(list(k=K, Q=quant[K], p=p))
#   
# }
# 
# 
# ProbMOM <- function(data, gamma, q, plot = FALSE, add = FALSE, 
#                     main = "Estimates of small exceedance probability", ...) {
#   
#   # Check input arguments
#   checkInput(data,gamma,gammapos=FALSE)
#   
#   if (length(q)>1) {
#     stop("q should be a numeric of length 1.")
#   }
#   
#   X <- as.numeric(sort(data))
#   n <- length(X)
#   prob <- numeric(n)
#   K <- 1:(n-1)
#   
#   # Estimator for tail probabilities
#   
#   H <- Hill(X)$gamma
#   a <- X[n-K]*H[K]*(1-pmin(gamma,0))
#   
#   prob[K] <- ((K+1)/(n+1)) * (1 + gamma[K]/a[K]*(q-X[n-K]))^(-1/gamma[K])
#   prob[prob<0 | prob>1] <- NA
#   
#   # plots if TRUE
#   plotfun(K, prob[K], type="l", xlab="k", ylab="1-F(x)", main=main, plot=plot, add=add, ...)
#   
#   # output list with values of k, corresponding return period estimates 
#   # and the considered large quantile q
#   
#   output(list(k=K, P=prob[K], q=q),plot=plot,add=add)
#   
# }

