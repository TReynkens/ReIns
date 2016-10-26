
########################################
# Compute estimate of 1-F_{Y|X}(y|x)
# censored indicates if an observation is uncensored
# y can be a vector
# x is a single number or a vector with the same length as y. In the latter case
# 1-F_{Y|X}(y|x) is computed for all n_y pairs (x,y) with n_y the length of y
crSurv <- function(x, y, Xtilde, Ytilde, censored, h, kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov")) {
  
  # Check input
  n <- length(Xtilde)
  
  if (!is.numeric(x)) {
    stop("x needs to be numeric.")
  }
  
  if (!is.numeric(y)) {
    stop("y needs to be numeric.")
  }
  
  if (!is.numeric(Xtilde)) {
    stop("Xtilde needs to be numeric.")
  }
  
  if (!is.numeric(Ytilde)) {
    stop("Ytilde needs to be numeric.")
  }
  
  if (length(Ytilde)!=n) {
    stop("Xtilde and Ytilde need to be of the same length.")
  }
  
  if (length(censored)!=n) {
    stop("Xtilde and censored need to be of the same length.")
  }
  
  # Give x always length(y)
  if (length(x)==1) {
    x <- rep(x, length(y))
    
  } else if (length(x)!=length(y)) {
    stop("x needs to have length 1 or the same length as y.")
  }
  
  # Set kernel to function with this name
  kernel <- match.arg(kernel)
  kernel <- .kernel_aux(kernel, h=1)
  
  # Convert censored to Delta
  Delta <- !as.logical(censored)
  

  # Create output vector
  surv <- numeric(length(y))
  

  # Sort data according to Ytilde values
  s <- sort(Ytilde, index.return=TRUE)
  Ytilde_sort <- s$x
  Xtilde_sort <- Xtilde[s$ix]
  Delta_sort <- Delta[s$ix]

  
  for (i in 1:length(y)) {
    
    # Nadaraya weights, only non-zero when uncensored
    # Sorted with Ytilde!
    weights_sort <- numeric(n)
    weights_sort[Delta_sort] <- kernel((x[i]-Xtilde_sort[Delta_sort])/h)
    weights_sort <- weights_sort / sum(weights_sort)

    # Select suitable indices
    ind <- which(Ytilde_sort<=y[i])
    
    # Compute cumulative sum of weights,
    # add 0 to avoid problems when ind is empty
    wsum <- c(0, cumsum(weights_sort[ind]))
    
    # Avoid numerical problems
    wsum[wsum==1] <- wsum[wsum==1] - 10^(-16)

    surv[i] <- prod((1-weights_sort[ind] / (1-wsum[-length(wsum)]))^Delta_sort[ind])
  }

  return(surv)
}

# Pareto QQ-plot adapted for censoring and regression
crParetoQQ <- function(x, Xtilde, Ytilde, censored, h, kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov"), 
                       plot = TRUE, add = FALSE, main = "Pareto QQ-plot", type = "p", ...) {
  

  if (length(x)!=1) {
    stop("x needs to have length 1.")
  }
  
  # Calculate theoretical and empirical quantiles
  
  surv <- crSurv(x=x, y=sort(Ytilde), Xtilde=Xtilde, Ytilde=Ytilde, censored=censored, 
                 h=h, kernel=kernel)
  # Remove too small values!
  pqq.the <- -log(surv[surv>=10^(-8)])
  pqq.emp <- log(sort(Ytilde)[surv>=10^(-8)])
  
  # plots if TRUE
  .plotfun(pqq.the, pqq.emp, type=type, xlab="-log(Survival probability)", ylab="log(Y)", 
           main=main, plot=plot, add=add, ...)
  
  # output list with theoretical quantiles pqq.the and empirical quantiles pqq.emp
  .output(list(pqq.the=pqq.the, pqq.emp=pqq.emp), plot=plot, add=add)
}


######################################################################
# Compute estimates of gamma using a Hill-type estimator
# censored indicates if an observation is censored
# x is a single number
crHill <- function(x, Xtilde, Ytilde, censored, h, kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov"), 
                   logk = FALSE, plot = FALSE, add = FALSE, main = "", ...) {

  # Check input
  n <- length(Xtilde)
  
  if (!is.numeric(x)) {
    stop("x needs to be numeric.")
  }
  
  if (!is.numeric(Xtilde)) {
    stop("Xtilde needs to be numeric.")
  }
  
  if (!is.numeric(Ytilde)) {
    stop("Ytilde needs to be numeric.")
  }
  
  if (length(Ytilde)!=n) {
    stop("Xtilde and Ytilde need to be of the same length.")
  }
  
  if (length(censored)!=n) {
    stop("Xtilde and censored need to be of the same length.")
  }
  
  
  if (length(x)!=1) {
    stop("x needs to have length 1.")
  }
  
  # Set kernel to function with this name
  kernel <- match.arg(kernel)
  kernel <- .kernel_aux(kernel, h=1)

  # Convert censored to Delta
  Delta <- !as.logical(censored)
  

  # Sort data according to Ytilde values
  s <- sort(Ytilde, index.return=TRUE)
  Ytilde_sort <- s$x
  Xtilde_sort <- Xtilde[s$ix]
  Delta_sort <- Delta[s$ix]
  
  # Nadaraya weights, only non-zero when uncensored
  # Sorted with Ytilde!
  weights_sort <- numeric(n)
  weights_sort[Delta_sort] <- kernel((x-Xtilde_sort[Delta_sort])/h)
  weights_sort <- weights_sort / sum(weights_sort)
  
  K <- 1:(n-1)
  H <- numeric(n)
  
  # Compute cumulative sum of weights,
  # add 0 to avoid problems when ind is empty
  wsum <- c(0, cumsum(weights_sort))
  
  # Avoid numerical problems
  wsum[wsum==1] <- wsum[wsum==1] - 10^(-16)
  
  for (k in 1:(n-1)) {

    prods <- numeric(k)
    
    # Compute products for i=1,...,k
    for (i in 1:k) {
      ind <- 1:(n-i)
      # Note that wsum[1]=0 !
      prods[i] <- prod((1 - weights_sort[ind] / (1-wsum[ind]))^Delta_sort[ind])
    }
    H[k] <- sum(prods * (log(Ytilde_sort[n-(1:k)+1])-log(Ytilde_sort[n-(1:k)]))) / prods[k]
  }
  
  # Plots
  if (logk) {
    .plotfun(log(K), H[K], type="l", xlab="log(k)", ylab="gamma", main=main, plot=plot, add=add, ...)
  } else {
    .plotfun(K, H[K], type="l", xlab="k", ylab="gamma", main=main, plot=plot, add=add, ...)
  }

  
  # output list with values of k and corresponding Hill estimates
  
  .output(list(k=K, gamma=H[K]), plot=plot, add=add)
}
