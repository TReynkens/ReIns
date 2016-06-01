
# Check input for const
.constCheck <- function(const) {
 
  if (!is.numeric(const)) stop("const should be numeric.")
  
  if (any(const<=0) | any(const>=1)) {
    stop("const should be a vector of numbers in (0,1).")
  }
  
  if (is.unsorted(const, strictly=TRUE)) {
    stop("const should be a strictly increasing vector.")
  }
  
}


# Check input for tsplice
.tspliceCheck <- function(tsplice) {
  
  if (!is.numeric(tsplice)) stop("tsplice should be numeric.")

  if (is.unsorted(tsplice, strictly=TRUE)) {
    stop("tsplice should be a strictly increasing vector.")
  }
  
}



#########################################################
# S3 classes


# Check if x is an integer
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
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
  if (!.is.wholenumber(M) | M<=0) stop("M must be a strictly positive integer.")
  
  # Check length of arguments
  if (length(p)!=M) stop("p must have length M.")
  if (length(shape)!=M) stop("shape must have length M.")
  if (length(theta)!=1) stop("theta must have length 1")
  
  eps <- .Machine$double.eps
  
  # Check if p is a vector of probabilities that sums to one.
  if (any(p<=eps)) stop("p must contain strictly positive numbers.")
  tol <- .Machine$double.eps^0.5
  if (abs(sum(p)-1)>tol) stop("p must have sum 1.")
  
  # Check if shape is a strictly increasing vector of positive integers.
  if (any(shape<=eps) | any(!.is.wholenumber(shape))) {
    stop("shape must contain strictly positive integers.")
  }  
  if (is.unsorted(shape, strictly=TRUE)) {
    stop("shape must be a strictly increasing vector.")
  }
  
  # Check if theta is strictly positive
  if (theta<=eps) {
    stop("theta must be a strictly positive number.")
  }
  
  # Make first list
  L <- list(p=p, shape=shape, theta=theta, M=M)
  
  # Optional argument M_initial
  if (!is.null(M_initial)) {
    
    if (!is.numeric(M_initial)) stop("M_initial must be numeric.")
    
    if (length(M_initial)!=1) stop("M_initial must have length 1.")
    
    if (!.is.wholenumber(M_initial) | M_initial<=0) {
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
    
    eps <- .Machine$double.eps
    if (any(sigma<=eps)) {
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
SpliceFit <- function(const, trunclower, t, type, MEfit, EVTfit, loglik = NULL, IC = NULL) {
  
  # Check input for const
  .constCheck(const)
  
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
  
  # Change tPa to TPa
  type[type=="tPa"] <- "TPa"
  
  if (length(type)!=l+1) stop("type should have one element more than const.")
  
  if (!type[1]=="ME") stop("The first argument of type is \"ME\".")
  
  if (!all(type[-1] %in% c("Pa", "GPD", "TPa"))) stop("Invalid type.")
  
  if (any(type[-1]=="GPD") & any(type[-1]!="GPD")) stop("GPD cannot be combined with other EVT distributions.")
  

  # Check MEfit and EVTfit
  if (class(MEfit)!="MEfit") stop("MEfit should be of class MEfit.")
  if (class(EVTfit)!="EVTfit") stop("EVTfit should be of class EVTfit.")
  
  if (length(EVTfit$gamma)!=l) stop("gamma should have the same length as const.")
  
  if (type[2]=="GPD" & !exists("sigma", where=EVTfit)) {
    stop("sigma should be included when GPD is used.")
  }

  NAind <- which(is.na(EVTfit$gamma))
  if (any(grepl("Pa", type[-c(1,NAind)]) & EVTfit$gamma[-NAind]<=0)) {
    stop("gamma should be strictly positive when Pareto distribution is used.")
  }
 

  # Check if correct type when truncated
  ind <- which(is.finite(EVTfit$endpoint))
  if (any(substring(type[ind+1],1,1)!="T")) {
    stop("Invalid type when endpoint is finite.")
  }
  
  # Check if truncated when indicated by type
  ind <- which(substring(type,1,1)=="T")
  if (any(!is.finite(EVTfit$endpoint[ind-1]) & !is.na(EVTfit$endpoint[ind-1]))) {
    stop("Infinite endpoint is not compatible with type \"TPa\".")
  }
  
  # Issue warning for NAs in gamma
  if (any(is.na(EVTfit$gamma))) {
    warning("One or more NAs are present in estimate(s) for gamma.")
  }
  
  # Make list
  L <- list(const=const, pi=pi, trunclower=trunclower, t=t, type=type, MEfit=MEfit, EVTfit=EVTfit)
  
  # Check log-likelihood and add to list
  if (!is.null(loglik)) {
    
    if (!is.numeric(loglik) | length(loglik)!=1) {
      stop("loglik should be a numeric of length 1.")
    }
    
    L$loglikelihood <- loglik
  }
  
  # Check IC and add to list
  if (!is.null(IC)) {
    
    if (!is.numeric(IC) | !(length(IC) %in% 1:3)) {
      stop("IC should be a numeric of length 1, 2 or 3.")
    }
    
    L$IC <- IC
  }
  
  # Return final structure of class Splicefit
  structure(L, class = "SpliceFit")
}

# Auxiliary function for summary.SpliceFit
.pasteVec <- function(s, name, vec, digits, paste_end="\n\n") {
  
  # Take no digits when vec is large
  digits <- ifelse(min(vec)<=1000, digits, 0)
  
  # Format output in right style (after rounding): add blank spaces, remove scientific notation and
  # remove leading whitespace, remove trailing zeros
  res <- format(round(vec,digits), big.mark=" ", scientific=FALSE, trim=TRUE,
                drop0trailing=TRUE)
  
  # paste_end is the last part of the string
  if(length(vec)==1) {
    paste0(s, name," = ",res, paste_end)
  } else {
    paste0(s, name," = (",paste(res, collapse=", "),")",paste_end)
  }
  
}



# summary method for SpliceFit class
# Needs to have same input name as generic one
summary.SpliceFit <- function(object, digits = 3, ...) {
  
  splicefit <- object
  mefit <- splicefit$MEfit
  evtfit <- splicefit$EVTfit
  
  # Length of dashes and stars
  l <- 20
  
  # Title
  s <- "\n"
  s <- paste0(s, paste0(rep("--",l),collapse=""), "\n")
  s <- paste0(s, "Summary of splicing fit", "\n")
  s <- paste0(s, paste0(rep("--",l),collapse=""), "\n\n")
  
  
  # General splicing part
  s <- .pasteVec(s, "const", splicefit$const, digits)
  s <- .pasteVec(s, "pi", splicefit$pi, digits)
  s <- .pasteVec(s, "t0", splicefit$trunclower, digits)
  s <- .pasteVec(s, "t", splicefit$t, digits)
  s <- paste0(s, "type = (",paste(splicefit$type, collapse=", "),")","\n\n")
  
  s <- paste0(s, paste0(rep("* ",l),collapse=""), "\n\n")
  
  
  # ME part
  s <- .pasteVec(s, "p", mefit$p, digits)
  s <- .pasteVec(s, "r", mefit$shape, digits)
  s <- .pasteVec(s, "theta", mefit$theta, digits)
  s <- .pasteVec(s, "M", mefit$M, digits)
  if (exists("M_initial", where=splicefit$MEfit)) {
    s <- .pasteVec(s, "M_initial", mefit$M_initial, digits)
  }
  
  s <- paste0(s, paste0(rep("* ",l),collapse=""), "\n\n")
  
  # EVT part
  s <- .pasteVec(s,"gamma", evtfit$gamma, digits)
  if (exists("sigma", where=splicefit$EVTfit)) {
    s <- .pasteVec(s,"sigma", evtfit$sigma, digits)
  }
  s <- .pasteVec(s,"endpoint", evtfit$endpoint, digits)
  
  cat(s)
}

# General tex function
tex <- function(x, ... ) UseMethod("tex", x)

# TeX method for SpliceFit class
tex.SpliceFit <- function(object, digits = 3, ...) {
  
  splicefit <- object
  mefit <- splicefit$MEfit
  evtfit <- splicefit$EVTfit
  
  s <- "\n"
  
  # General splicing part
  
  # Print all but last element of pi
  if(length(splicefit$pi)>2) {
    for(i in 1:(length(splicefit$pi)-1)) {
      s <- .pasteVec(s, paste0("\\pi_",i), splicefit$pi[i], digits, ", ")
    }
    # Remove last ", " and add newline (twice)
    s <- substr(s,1,nchar(s)-2)
    s <- paste0(s, "\n\n")
    
  } else {
    
    s <- .pasteVec(s, "\\pi", splicefit$pi[1], digits)
  }
  
  s <- .pasteVec(s, "t^l", splicefit$trunclower, digits)
  
  # Print all splicing points
  if(length(splicefit$t)>1) {
    
    for(i in 1:length(splicefit$t)) {
      s <- .pasteVec(s, paste0("t_",i), splicefit$t[i], digits, ", ")
    }
    # Remove last ", " and add newline (twice)
    s <- substr(s,1,nchar(s)-2)
    s <- paste0(s, "\n\n")
    
  } else {
    s <- .pasteVec(s, "t", splicefit$t[1], digits)
  }
  
  # Print splicing type
  s <- paste0(s, "type = (",paste(splicefit$type, collapse=", "),")","\n\n")
  
  
  # ME part
  # Print alpha
  s <- .pasteVec(s, "\\alpha", mefit$p, digits)
  # Print shape
  s <- .pasteVec(s, "r", mefit$shape, digits)
  # Print scale
  s <- .pasteVec(s, "\\theta", mefit$theta, digits)
  
  # EVT part
  
  # Print gamma
  if(length(evtfit$gamma)>1) {
    
    for(i in 1:length(evtfit$gamma)) {
      s <- .pasteVec(s, paste0("\\gamma_",i), evtfit$gamma[i], digits, ", ")
    }
    # Remove last ", " and add newline (twice)
    s <- substr(s,1,nchar(s)-2)
    s <- paste0(s, "\n\n")
    
  } else {
    s <- .pasteVec(s, "\\gamma", evtfit$gamma, digits)
  }
  
  # Print sigma if exists
  if (exists("sigma", where=splicefit$EVTfit)) {
    if(length(evtfit$sigma)>1) {
      
      for(i in 1:length(evtfit$sigma)) {
        s <- .pasteVec(s, paste0("\\sigma_",i), evtfit$sigma[i], digits, ", ")
      }
      # Remove last ", " and add newline (twice)
      s <- substr(s,1,nchar(s)-2)
      s <- paste0(s, "\n\n")
      
    } else {
      s <- .pasteVec(s, "\\sigma", evtfit$sigma, digits)
    }
  }
  
  # Print endpoint, special case when infinite
  end <- evtfit$endpoint[length(evtfit$endpoint)]
  
  if (is.finite(end)) {
    s <- .pasteVec(s,"T", end, digits)
  } else {
    s <- paste0(s, "T = +\\infty","\n")
  }
  
  
  cat(s)
}



# Only keep necessary parts from .MEtune output
.MEoutput <- function(fit_tune) {
  
  # Obtain best model
  MEfit_old <- fit_tune$best_model
  
  # Make MEfit object
  # Use alpha not beta for future calculations!
  MEfit <- MEfit(p=MEfit_old$alpha, shape=MEfit_old$shape, theta=MEfit_old$theta, 
                 M=MEfit_old$M, M_initial=MEfit_old$M_initial)

  return(MEfit)
}


###############################################################################

# Fit splicing of mixed Erlang and (truncated) Pareto
SpliceFitPareto <- function(X, const = NULL, tsplice = NULL, M = 3, s = 1:10, trunclower = 0, truncupper = Inf, EVTtruncation = FALSE, 
                            ncores = NULL, criterium = c("BIC","AIC"), reduceM = TRUE, 
                            eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
 
  ##
  # Check input
  
  # Check if X is numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  n <- length(X)
  
  # Check input for const if not NULL
  if (!is.null(const)) {

    .constCheck(const)
    
    l <- length(const)
  }
  
  # Check input for tsplice if not NULL
  if (!is.null(tsplice)) {
    
    .tspliceCheck(tsplice)
    
    l <- length(tsplice)
  }

  if (is.null(const) & is.null(tsplice)) {
    stop("const and tsplice cannot be both NULL.")
  }
  
  if (!is.null(const) & !is.null(tsplice)) {
    warning("Both const and tsplice are not NULL, the input for tsplice will be used.")
    
    const <- NULL
  }
  
  
  # Check input for EVTtruncation
  if (!is.logical(EVTtruncation)) {
    stop("EVTtruncation should be a logical of length 1.")
  }
  
  if (length(EVTtruncation)!=1 & length(EVTtruncation)>0) {
    # Only use last element if length>1 (compatibility with older versions)
    EVTtruncation <- EVTtruncation[length(EVTtruncation)]
    warning("EVTtruncation has more than one element, only the last one is used.")
    
  } else if(length(EVTtruncation)<=0) {
    # Stop if length 0 (or less)
    stop("EVTtruncation should be a logical of length 1.")
  }
  
  # Check input for trunclower
  if (length(trunclower)!=1) {
    stop("trunclower should have length 1.")
  }
  
  if (!is.numeric(trunclower)) {
    stop("trunclower should be numeric.")
  }
  
  if (trunclower>min(X)) {
    stop("trunclower should be smaller than all data points.")
  }
  
  if (trunclower<0) {
    stop("trunclower cannot be strictly negative.")
  }
  
  if (!is.finite(trunclower) ) {
    stop("trunclower has to be finite.")
  }
  
  
  # Check input for truncupper
  if (length(truncupper)!=1 ) {
    stop("truncupper should have length 1.")
  }
  
  if (!is.numeric(truncupper)) {
    stop("trunclower should be numeric.")
  }
  
  if (truncupper<max(X)) {
    stop("truncupper should be larger than all data points.")
  }

  if (truncupper<=trunclower) {
    stop("truncupper should be strictly larger than trunclower.")
  }
  
  if (!EVTtruncation & is.finite(truncupper)) {
    # Set EVTtruncation to true when truncupper is finite
    EVTtruncation <- TRUE
    warning("truncupper is finite, EVTtruncation is set to TRUE.")
  }
  
  
  
  # Match input argument for criterium
  criterium <- match.arg(criterium)
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  
  ##
  # Determine values of splicing points
    
  Xsort <- sort(X)
  
  if (!is.null(const)) {
    # Number of points larger than splicing points
    k_init <- ceiling((1-const)*n)
    
    # Splicing points
    tvec <- numeric(l)
    tvec <- Xsort[n-k_init]
    
    # Update const
    const <- 1-k_init/n

  } else {
    tvec <- tsplice
    
    # const is number of observations smaller than or equal to tvec
    f <- function(t, x) sum(x<=t)
    const <- sapply(tvec, f, Xsort)/n

    # Make sure k_init is of integer type
    k_init <- as.integer((1-const)*n)
  }

  # Number of points in splicing part larger than splicing point
  kvec <- k_init
  
  if (l>1) {
    for(i in (l-1):1) {
      kvec[i] <- kvec[i] - sum(kvec[(i+1):l])
    } 
  }
  

  # Problem when first splicing point smaller than trunclower
  if (any(trunclower>=tvec[1])) {
    stop("trunclower should be strictly smaller than the first splicing point.")
  }
  
  t1 <- tvec[1]
   
  ##
  # Mixing Erlang part
  MEind <- (X<=t1) 
  
  # Upper truncated at threshold t
  fit_tune <- .MEtune(lower=X[MEind], upper=X[MEind], trunclower=trunclower, truncupper=t1,
                      M=M, s=s, nCores=ncores, criterium=criterium, reduceM=reduceM, 
                      eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  # Output as MEfit object
  MEfit <- .MEoutput(fit_tune)
  
  ##
  # EVT part
  EVTfit <- list()
  
  type <- character(l)
  
  EVTfit$gamma <- numeric(l)
  EVTfit$endpoint <- rep(Inf, l)

  # Splicing parts
  for (i in 1:l) {

    if (i==l) {
      
      if (EVTtruncation) {
        
        if (is.finite(truncupper)) {
          # Last Pareto distribution is truncated, known truncation point
          EVTfit$gamma[i] <- .trHillinternal(X[X<=truncupper], tvec[i], truncupper)
          EVTfit$endpoint[i] <- truncupper
          
        } else {
          # Last Pareto distribution is truncated, estimate truncation point
          res <- trHill(X)
          resDT <- trDT(X, gamma=res$gamma)
          resEndpoint <- trEndpoint(X, gamma=res$gamma, DT=resDT$DT)
          EVTfit$gamma[i] <- res$gamma[res$k==kvec[i]]
          EVTfit$endpoint[i] <- resEndpoint$Tk[resEndpoint$k==kvec[i]]
        }

        type[i] <- "TPa"
        
      } else {
        # Last Pareto distribution is not truncated
        
        EVTfit$gamma[i] <- .Hillinternal(X, tvec[i])
        type[i] <- "Pa"
      }
  
    } else {
      
      # Endpoint is next splicing point
      tt <- tvec[i+1]
      
      # Truncated Pareto distribution with truncation point tt
      EVTfit$gamma[i] <- .trHillinternal(X[X<=tt], tvec[i], tt)
      EVTfit$endpoint[i] <- tt
      type[i] <- "TPa"
    }
    
  }
  
  # Convert to object of class EVTfit
  EVTfit <- structure(EVTfit, class="EVTfit")
  
  
  ##

  # Make SpliceFit object
  sf <- SpliceFit(const=const, trunclower=trunclower, t=tvec, type=c("ME",type), MEfit=MEfit, EVTfit=EVTfit)
  
  # Compute log-likelihood
  loglik <- sum(log(dSplice(X, sf)))

  # Compute ICs, parameters: alpha (-1 since sums to 1), r, theta, gamma's, last endpoint (if finite), pi's
  df <- 2*MEfit$M-1 + 1 + length(EVTfit$gamma) + is.finite(EVTfit$endpoint[length(EVTfit$endpoint)]) + length(const)
  aic <- -2*loglik + 2*df
  bic <- -2*loglik + log(n)*df
  ic <- c(aic, bic)
  names(ic) <- c("AIC", "BIC")
  
  # Add to SpliceFit object
  sf$loglik <- loglik
  sf$IC <- ic
  
  ##
  # Return SpliceFit object
  return( sf )
}
# Include for compatibility with old versions
SpliceFitHill <- SpliceFitPareto


# Fast Hill estimator for fixed threshold
.Hillinternal <- function(x, threshold) {
  
  ind <- (x>threshold)
  
  return(sum(log(x[ind]/threshold))/sum(ind))
}

# Fast Hill estimator for fixed threshold and known truncation point 
.trHillinternal <- function(x, threshold, truncupper) {
  
  H <- .Hillinternal(x, threshold)
  
  beta <- truncupper/threshold
  f <- function(gamma) (gamma-H-log(beta)/(beta^(1/gamma)-1))^2
  
  return(nlm(f,H)$estimate)
}


##################################

# Fit splicing of ME and Pareto distribution to interval censored data
SpliceFitciPareto <- function(L, U, censored, tsplice, M = 3, s = 1:10, trunclower = 0, truncupper = Inf, ncores = NULL, 
                              criterium = c("BIC","AIC"), reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {
  
  warning("This function has not yet been implemented.")
}


##################################

# Fit splicing of mixed Erlang and GPD (POT)
# No truncation implemented, so only use with one GPD part!
SpliceFitGPD <- function(X, const = NULL, tsplice = NULL, M = 3, s = 1:10, trunclower = 0, ncores = NULL, 
                         criterium = c("BIC","AIC"), reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf) {

  # Check if X is numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  n <- length(X)
  
  
  # Check input for const if not NULL
  if (!is.null(const)) {
    
    .constCheck(const)
    
    l <- length(const)
  }
  
  # Check input for tsplice if not NULL
  if (!is.null(tsplice)) {
    
    .tspliceCheck(tsplice)
    
    l <- length(tsplice)
  }
  
  if (is.null(const) & is.null(tsplice)) {
    stop("const and tsplice cannot be both NULL.")
  }
  
  if (!is.null(const) & !is.null(tsplice)) {
    warning("Both const and tsplice are not NULL, the input for tsplice will be used.")
    
    const <- NULL
  }
  
  # Only implemented for a single splicing point
  if (l>=2) {
    stop("const (or tsplice) should be a single number.")
  }
  
  
  # Check input for trunclower
  if (!is.numeric(trunclower)) {
    stop("trunclower should be numeric.")
  }
  
  if (any(trunclower>min(X))) {
    stop("trunclower should be strictly smaller than all data points.")
  }
  
  if (any(trunclower<0)) {
    stop("trunclower cannot be strictly negative.")
  }
  
  
  
  # Check input for ncores
  if (is.null(ncores)) ncores <- max(detectCores()-1, 1)
  if (is.na(ncores)) ncores <- 1
  
  # Match input argument for criterium
  criterium <- match.arg(criterium)

  
  ##
  # Determine values of splicing points
  
  Xsort <- sort(X)
  
  if (!is.null(const)) {
    # Number of points larger than splicing points
    k_init <- ceiling((1-const)*n)
    
    # Splicing points
    tvec <- numeric(l)
    tvec <- Xsort[n-k_init]
    
    # Update const
    const <- 1-k_init/n
    
  } else {
    tvec <- tsplice
    
    # const is number of observations smaller than or equal to tvec
    f <- function(t, x) sum(x<=t)
    const <- sapply(tvec, f, Xsort)/n
    
  }
  
  
  # Problem when first splicing point smaller than trunclower
  if (any(trunclower>=tvec[1])) {
    stop("trunclower should be strictly smaller than the first splicing point.")
  }
  
  t1 <- tvec[1]
  
  ##
  # Mixing Erlang part
  MEind <- (X<=t1) 
  
  # Upper truncated at threshold t
  fit_tune <- .MEtune(lower=X[MEind], upper=X[MEind], trunclower=trunclower, truncupper=t1,
                     M=M, s=s, nCores=ncores, criterium=criterium, reduceM=reduceM, 
                     eps=eps, beta_tol=beta_tol, maxiter=maxiter)
  # Output as MEfit object
  MEfit <- .MEoutput(fit_tune)
  
  ##
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
  
  # Make SpliceFit object
  sf <- SpliceFit(const=const, trunclower=trunclower, t=tvec, type=c("ME",type), MEfit=MEfit, EVTfit=EVTfit)
  
  # Compute log-likelihood
  loglik <- sum(log(dSplice(X, sf)))
  
  # Compute ICs, parameters: alpha (-1 since sums to 1), r, theta, gamma's and sigma's, last endpoint (if finite), pi's
  df <- 2*MEfit$M-1 + 1 + 2*length(EVTfit$gamma) + is.finite(EVTfit$endpoint[length(EVTfit$endpoint)]) + length(const)
  aic <- -2*loglik + 2*df
  bic <- -2*loglik + log(n)*df
  ic <- c(aic, bic)
  names(ic) <- c("AIC", "BIC")
  
  # Add to SpliceFit object
  sf$loglik <- loglik
  sf$IC <- ic
  
  
  ##
  # Return SpliceFit object
  return( sf )
}



########################################################################

# Splicing PDF
dSplice <- function(x, splicefit, log = FALSE) {
  
  # Check input
  if (class(splicefit)!="SpliceFit") stop("splicefit should be of class SpliceFit.")
  
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  
  d <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  if (any(is.na(EVTfit$gamma))) {
    stop("At least one of the estimates for gamma is not available.")
  }  
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<tvec[1])
  
  # Case x<=t
  d[ind] <- const[1] * .ME_density(x[ind], shape = MEfit$shape, alpha = MEfit$p, 
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
      
    } else if (type[i+1] %in% c("Pa","TPa")) {
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
  
  # Check input
  if (class(splicefit)!="SpliceFit") stop("splicefit should be of class SpliceFit.")
  
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  
  p <- numeric(length(x))
  
  MEfit <- splicefit$MEfit
  EVTfit <- splicefit$EVTfit
  
  if (any(is.na(EVTfit$gamma))) {
    stop("At least one of the estimates for gamma is not available.")
  }  
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)
  
  type <- splicefit$type
  
  ind <- (x<tvec[1])
  
  # Case x<t
  p[ind] <- const[1] * .ME_cdf(x[ind], shape = MEfit$shape, alpha = MEfit$p, 
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
    
    } else if (type[i+1] %in% c("Pa","TPa")) {
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
  if (class(splicefit)!="SpliceFit") stop("splicefit should be of class SpliceFit.")

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
  
  if (any(is.na(EVTfit$gamma))) {
    stop("At least one of the estimates for gamma is not available.")
  }  
  
  tvec <- splicefit$t
  trunclower <- splicefit$trunclower
  
  const <- splicefit$const
  l <- length(const)

  ind <- (p<const[1])
  
  if (any(ind)) {
    # Quantiles of ME part
    q[ind] <- .ME_VaR(p[ind]/const[1], shape = MEfit$shape, alpha = MEfit$p, 
                     theta = MEfit$theta, trunclower=trunclower, truncupper=tvec[1], 
                     interval=c(trunclower,tvec[1])) 
  }
  
  
  # Quantiles of EVT part
  
  for (i in 1:l) {
    
    # Next value for const (1 for last part)
    cconst <- ifelse(i==l, 1, const[i+1])
    
    # Index for all probabilities in i-th EVTpart
    ind <- p>=const[i] & p<cconst
    
    # Next splicing point (Inf for last part)
    tt <- ifelse(i==l, Inf, tvec[i+1])
    e <- min(EVTfit$endpoint[i], tt)
    
    # i+1 since type[1]="ME"
    if (splicefit$type[i+1]=="GPD") {
      q[ind] <- qtgpd((p[ind]-const[i])/(cconst-const[i]), mu=tvec[i], gamma=EVTfit$gamma[i], sigma=EVTfit$sigma[i], endpoint=e)
      
    } else if (splicefit$type[i+1] %in% c("Pa","TPa")) {
      q[ind] <- qtpareto((p[ind]-const[i])/(cconst-const[i]), shape=1/EVTfit$gamma[i], scale=tvec[i], endpoint=e)
      
    } else {
      stop("Invalid type.")
    }
  }
  
  # Special case
  q[p==1] <- EVTfit$endpoint[l]
  
  return(q)
}

# Random numbers from spliced distribution
rSplice <- function(n, splicefit) {
  
  return(qSplice(runif(n), splicefit=splicefit))
}

###########################################################################



# Plot of fitted survival function and ECDF estimator + bounds
SpliceECDF <- function(x, X, splicefit, alpha = 0.05, ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
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
SpliceTB <- function(x = sort(L), L, U = L, censored, splicefit, alpha = 0.05, ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Plot fitted survival function
  plot(x, 1-pSplice(x, splicefit=splicefit), type="l", xlab="x", ylab="1-F(x)", ...)
  
  
  # Add Turnbull survival function
  tb <- Turnbull(x, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                 truncupper=max(splicefit$EVTfit$endpoint), conf.type="plain", conf.int=1-alpha)
  lines(x, tb$surv, col="red")
  
  # Add confidence intervals
  fit <- tb$fit
  lines(fit$time, fit$lower, col = "blue", lty = 2)
  lines(fit$time, fit$upper, col = "blue", lty = 2)
  legend("topright", c("Fitted survival function","Turnbull estimator","95% confidence intervals"),
         lty=c(1,1,2), col=c("black","red","blue"))
}


# Probability - probability plot with ECDF
SplicePP <- function(X, splicefit, x = sort(X), log = FALSE, ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # ECDF estimator
  fit  <- ecdf(X)
  est <- 1-fit(x)
  
  # Plot fitted survival function vs. empirical survival function or use minus log-versions
  if (log) {
    ind <- est>0
    plot(-log(est[ind]), -log(1-pSplice(x[ind],splicefit=splicefit)), type="p",
         xlab="-log(Empirical survival probability)",
         ylab="-log(Fitted survival probability)", ...)
  } else {
    plot(est, 1-pSplice(x,splicefit=splicefit), type="p", xlab="Empirical survival probability",
         ylab="Fitted survival probability", ...)
  }
  abline(a=0, b=1)
}

# Probability - probability plot with Turnbull estimator
SplicePP_TB <- function(L, U = L, censored, splicefit, x = NULL, log = FALSE, ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.null(x) & !is.numeric(x)) stop("x should be a numeric vector or NULL.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                  truncupper=max(splicefit$EVTfit$endpoint))
    if (is.null(x)) {
      # # Use unique points of Turnbull intervals since Turnbull estimator is exact there
      # x <- unique(as.numeric(SurvTB$fit$intmap))
      # Use empirical quantiles
      n <- length(L)
      p <- (1:n) / (n+1)
      x <- SurvTB$fquant(p)
    }
    
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                     truncupper=max(splicefit$EVTfit$endpoint))
    
    if (is.null(x)) {
      # # Use knots (jump points)
      # x <- unique(SurvTB$fit$knots)
      # Use empirical quantiles
      n <- length(L)
      p <- (1:n) / (n+1)
      x <- SurvTB$fquant(p)
    }
  }

  # Survival function in x
  surv <- SurvTB$f(x)

  # Plot fitted survival function vs. Turnbull survival function or use minus log-versions
  if (log) {
    ind <- surv>0
    plot(-log(surv[ind]), -log(1-pSplice(x[ind],splicefit=splicefit)), type="p", 
         xlab="-log(Turnbull survival probability)",
         ylab="-log(Fitted survival probability)", ...)
  } else {
    plot(surv, 1-pSplice(x, splicefit=splicefit), type="p",
         xlab="Turnbull survival probability", ylab="Fitted survival probability", ...)
  }
  abline(a=0,b=1)
}



# Quantile - quantile plot with empirical quantiles
SpliceQQ <- function(X, splicefit, p = NULL, plot = TRUE, main = "Splicing QQ-plot", ...) {
  
  # Check input arguments
  .checkInput(X, pos=FALSE)
  
  X <- as.numeric(sort(X))
  n <- length(X)
  
  # Probabilities
  if (is.null(p)) {
    p <- (1:n)/(n+1)
    
  } else {
    
    if (!is.numeric(p)) {
      stop("p should be numeric.")
    }
    
    if(is.infinite(max(splicefit$EVTfit$endpoint)) & any(p==1)) {
      stop("All elements of p should be strictly smaller than 1 since the splicing
           distribution has an infinite endpoint.")
    }
      
    if (any(p<0 | p>1)) {
      stop("All elements of p should be between 0 and 1.")
    }
    
    if (length(p)!=n) {
      stop("p should have the same length as x.")
    }
    
    p <- sort(p)
  }
  
  
  # Empirical quantiles
  sqq.emp <- sort(X)
  
  # Quantiles of fitted distribution
  sqq.the <- qSplice(p=p, splicefit=splicefit)
  
  .plotfun(sqq.the, sqq.emp, type="p", xlab="Quantiles of splicing fit", ylab="Empirical quantiles", 
           main=main, plot=plot, add=FALSE, ...)
  
  # Add 45 degree line
  if (plot) abline(0,1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(sqq.the=sqq.the, sqq.emp=sqq.emp), plot=plot, add=FALSE)
}


# QQ-plot with Turnbull estimator
SpliceQQ_TB <- function(L, U = L, censored, splicefit, plot = TRUE, main = "Splicing QQ-plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }

  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE) & !all(censored==rep(0, length(L)))) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                  truncupper=max(splicefit$EVTfit$endpoint))

    # # Use unique points of Turnbull intervals since Turnbull estimator is exact there
    # x <- unique(as.numeric(SurvTB$fit$intmap))
    # # Corresponding function values
    # s <- SurvTB$f(x)
    # p <- 1-s

  } else {

    # Special warning when no censoring
    if (all(censored==rep(0, length(L)))) {
      warning("Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    } else {
      warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
              call.=FALSE)
    }
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                                 truncupper=max(splicefit$EVTfit$endpoint))
    
    # # Use knots (jump points)
    # x <- knots(SurvTB$f)
    # # Corresponding function values
    # s <- eval(expression(y), envir = environment(SurvTB$f))
    # p <- 1-s
  }

  
  # Probabilities
  n <- length(L)
  p <- (1:n) / (n+1)
  # Use empirical quantiles
  x <- SurvTB$fquant(p)

  # Remove 1 if infinite endpoint
  if(is.infinite(max(splicefit$EVTfit$endpoint))) {
    sqq.emp <- x[p<1-10^(-5)]
    p <- p[p<1-10^(-5)]
  } else {
    sqq.emp <- x
  }
  
  # Quantiles of fitted distribution
  sqq.the <- qSplice(p=p, splicefit=splicefit)
  
  .plotfun(sqq.the, sqq.emp, type="p", xlab="Quantiles of splicing fit", ylab="Empirical quantiles", 
             main=main, plot=plot, add=FALSE, ...)
  
  # Add 45 degree line
  if (plot) abline(0,1)
  
  # output list with theoretical quantiles sqq.the and empirical quantiles sqq.emp
  .output(list(sqq.the=sqq.the, sqq.emp=sqq.emp), plot=plot, add=FALSE)
}





# Log-log plot with empirical survival function and fitted survival function
SpliceLL <- function(x = sort(X), X, splicefit, ...) {
  
  # Check if X and x are numeric
  if (!is.numeric(X)) stop("X should be a numeric vector.")
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # ECDF estimator
  fit  <- ecdf(X)
  X <- sort(X)
  
  # Plot log of empirical survival function vs. sorted values
  plot(log(X), log(1-fit(X)), ylab="log(empirical survival probability)", xlab="log(X)", type="p", ...)
  # Add log of fitted survival function
  lines(log(x), log(1-pSplice(x, splicefit=splicefit)))
}


# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB <- function(x = sort(L), L, U = L, censored, splicefit, ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check if x is numeric
  if (!is.numeric(x)) stop("x should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  
  Zs <- sort(L)
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull2(Zs, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                         truncupper=max(splicefit$EVTfit$endpoint))
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    SurvTB <- Turnbull(Zs, L=L, R=U, censored=censored, trunclower=splicefit$trunclower,
                       truncupper=max(splicefit$EVTfit$endpoint))
  }

  # Plot log of Turnbull survival function vs. sorted values
  plot(log(Zs), log(SurvTB$surv), ylab="log(Turnbull survival probability)", xlab="log(X)", type="p", ...)
  # Add log of fitted survival function
  lines(log(x), log(1-pSplice(x, splicefit=splicefit)))
}


# Mean excess plot using Turnbull estimator
MeanExcess_TB <- function(L, U = L, censored, trunclower = 0, truncupper = Inf, 
                          plot = TRUE, k = TRUE, main = "Mean excess plot", ...) {
  
  # Check if L and U are numeric
  if (!is.numeric(L)) stop("L should be a numeric vector.")
  if (!is.numeric(U)) stop("U should be a numeric vector.")
  
  # Check lengths
  if (length(L)!=length(U)) stop("L and U should have the same length.")
  
  if (length(censored)==1) censored <- rep(censored, length(L))
  
  if (length(censored)!=length(L)) {
    stop("censored should have length 1 or the same length as L and U.")
  }
  
  # Turnbull survival function
  if (requireNamespace("interval", quietly = TRUE)) {
    SurvTB <- .Turnbull_internal2(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$xall
    y <- SurvTB$survall
    m <- length(x)
    
  } else {
    warning("Package \"interval\" is not available, Turnbull survival function from the \"survival\" package is used.", 
            call.=FALSE)
    
    SurvTB <- .Turnbull_internal(L=L, R=U, censored=censored, trunclower=trunclower, truncupper=truncupper)
    
    x <- SurvTB$fit$time
    y <- SurvTB$fit$surv
    m <- length(x)
  }
  
  
  # Values to evaluate mean excess function in
  n <- length(L)
  p <- (1:n) / (n+1)
  ic <- SurvTB$fquant(p)
  
  # Keep unique elements
  eps <- sqrt(.Machine$double.eps)
  ic[which(diff(ic)<=eps)+1] <- NA
  ic <- ic[!is.na(ic)]
  
  
  n <- length(ic)
  K <- 1:(n-1)
  me <- numeric(n)
  me2 <- numeric(n)
  
  # Compute slopes and intercepts of function in intervals
  as <- diff(y) / diff(x)
  as[diff(x)<=eps] <- 0
  bs <- y[1:(m-1)] - as * x[1:(m-1)]
  # Remove influence of jumps
  bs[is.infinite(as)] <- 0
  as[is.infinite(as)] <- 0
  
  for (i in K) {
    
    # Threshold
    r <- ic[n-i]
    
    if (r>x[m]) {
      me[i] <- 0
      
    } else {
      
      # Find correct interval (right-closed)
      ind_r <- length(x) - findInterval(-r, -rev(x))
      #ind_r = findInterval(r, x)
      
      if (r<x[1]) r <- x[1]
      
      # Slope and intercept of function in interval containing r
      if(ind_r==0) {
        # Special case if ind_r <- 0
        a <- 0
        b <- 1
        
      } else {
        a <- as[ind_r]
        b <- bs[ind_r]
      }
      
      # Use the fact that the Turnbull survival function is piecewise linear with jumps
      # First, contribution of interval containing r, then contributions of all intervals after this interval
      me[i] <- a * (x[ind_r+1]^2-r^2) / 2 + b *(x[ind_r+1]-r)  + 
        sum(as[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]^2) / 2 + bs[(ind_r+1):(m-1)] * diff(x[(ind_r+1):m]))
    }
    
  }
  
  # Mean-excess is int_r^Inf (1-F(z)) dz / (1-F(r))
  me[K] <- me[K] / SurvTB$f(ic[n-K])
  
  # Plot estimates
  if (plot) {
    
    if (k) {
      .plotfun(K, me[K], type="p", xlab="k", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    } else {
      .plotfun(ic[n-K], me[K], xlab=bquote(X["n-k,n"]), type="p", ylab=bquote(e["k,n"]), main=main, plot=TRUE, add=FALSE, ...)
    }
    
  }
  
  # Output list with values of k, order statistics X_n-k,n 
  # and mean excess scores e_k,n
  .output(list(k=K, X=ic[n-K], e=me[K]), plot=plot, add=FALSE)
}

