
###################################
#Several auxiliary functions
###################################

# Auxiliary function to handle plots for all estimators
.plotfun <- function(x, y, type, xlab, ylab, main, plot, add, ... ){
  
  # Handle the case that both plot and add are TRUE
  if (plot & add) {
    
    # returns FALSE if plot.new has not been called (which generated a warning in par(new=TRUE))
    # otherwise we return TRUE (no warning in par(new=TRUE))
    called <- tryCatch( {par(new=TRUE); TRUE}, warning = function(x) FALSE)
    
    if (!called) {
      # plot.new has not been called, create new plot
      warning("Both \'plot\' and \'add\' are TRUE, \'add\' is set to FALSE: 
              the results are plotted in a separate plot.")
      add <- FALSE
    } else {
      # plot.new has been called, add to this plot
      warning("Both \'plot\' and \'add\' are TRUE, \'plot\' is set to FALSE: 
              the results are added to an existing plot.")
      plot <- FALSE
    }
 
  }
  
  if (plot | add) {
    if ( plot ) { 
      # plot estimates
      plot(x, y, type=type, xlab=xlab, ylab=ylab, main=main, ...)
    }  else { 
      # add estimates to existing plot
      lines(x, y, ...)
    }
  }
}


####################################################

# Check input for estimators; gamma, scale and DT are optional arguments
# pos indicates if it should be checked that all elements of data
#     are strictly positive
# gammapos indicates if it should be checked that all elements of gamma
#     are strictly positive
# scalepos indicates if it should be checked that all elements of scale
#     are strictly positive
# DTpos indicates if it should be checked that all elements of DT
#     are positive
# r is a parameter for truncation
.checkInput <- function (data, gamma, scale, DT, pos = TRUE, gammapos = TRUE,
                        scalepos = TRUE, DTpos = TRUE, r = 1){
  
  ##################
  # data
  
  # check if data is numeric
  if (!is.numeric(data)) {
    stop("data should be a numeric vector.")
  }
    
  n <- length(data)
  
  # check if we have at least two data points
  if (n==1) {
    stop("We need at least two data points.")
  }
  
  if (pos) {
    # check if all elements of data are strictly positive
    if (min(data)<=0) {
      stop("data can only contain strictly positive values.")
    }
  }
  
  ##################
  # gamma

  # Check conditions for gamma if gamma is given
  if (!missing(gamma)) {
    
    if (!is.numeric(gamma)) {
      stop("gamma should be a numeric vector.")
    }
    
    if (gammapos & min(gamma,na.rm=TRUE)<=0) {
        stop("gamma can only contain strictly positive values.")
    }
    
    if (!(length(gamma) %in% c(n-2,n-1,n) + (r-1) )) {
      stop(paste0("gamma should have length ",n-2,", ",n-1," or ",n,"."))
    }
  } 
  
  
  ##################
  # scale
  
  # Check conditions for scale if scale is given
  if (!missing(scale)) {
    
    if (!is.numeric(scale)) {
      stop("scale should be a numeric vector.")
    }
    
    if (scalepos & min(scale,na.rm=TRUE)<=0) {
      stop("scale can only contain strictly positive values.")
    }
    
    if (!(length(scale) %in% c(n-2,n-1,n))) {
      stop(paste("scale should have length",n-2,",",n-1,"or",n))
    }
  }  
  
  if(!missing(gamma) & !missing(scale)) {
    if(length(gamma) != length(scale)) {
      stop("gamma and scale should have the same length.")
    }
  }
  
  
  ##################
  # DT
  
  # Check conditions for DT if DT is given
  if (!missing(DT)) {
    
    # check if DT is numeric
    if (!is.numeric(DT)) {
      stop("DT should be a numeric vector.")
    }
    
    if (DTpos) {
      # check if all elements of DT are positive
      if (min(DT,na.rm=TRUE)<0) {
        stop("DT can only contain positive values.")
      }
    }
    
    #Check length
    if(length(DT) !=1 & length(DT)!=(n-r)) {
      stop(paste("DT should have length 1 or length", n-r))
    }
  }
}

# Check if p is a vector of probabilities of length l
.checkProb <- function(p, l = 1) {
  
  # Check length
  if (length(p)!=l) {
    stop(paste0("p should be a numeric of length ",l, "."))
  }
  
  # Check if valid probability
  if (is.numeric(p)) {
    if(p<0 | p>1) {
      stop("p should be between 0 and 1.")
    }
  } else {
    stop("p should be numeric.")
  }
}


# Check if censored is a logical of the right length
# n is the length of the data
.checkCensored <- function(censored, n) {
  
  if (length(censored)!=1) {
    if (n != length(censored)) {
      stop("data and censored should have the same length.")
    }
  } else {
    censored = rep(censored,n)
  }

  # Check if logical vector or a vector with only 0 and 1
  if (!is.logical(censored)) {
    if (!all(censored ==1 | censored==0)) {
      stop("censored should be a logical vector.")
    } else {
      censored <- as.logical(censored)
    }
    
  }
  
  
  # Check if all are censored
  if (all(censored)) {
    stop("Not all data points can be censored.")
  }
  
  return(censored)
}


####################################################

# If a plot is made or results are added to a plot, x is not printed
# but x can still be assigned!
.output <- function(x, plot, add) {
  
  if (plot || add) {
    return(invisible(x))
  } else {
    return(x)
  }
}


####################################################

