# Call C++ functions for fitting the ME-Pa splicing model.
# Some input management is done here.
# This file should not be confused with "RcppExports.R"!

.spliceEM_fit_raw_cpp <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf,
                                  pi, theta, shape, beta, gamma,
                                  eps = 10 ^ (-3), beta_tol = 10 ^ (-5), maxiter = Inf) {

  # Separate uncensored and censored observations
  uncensored <- !censored

  # Indices for the five different cases
  # Some can be numeric(0)!

  # t^l <= x_i=l_i=u_i <= t <T
  ind1 <- which(uncensored & (upper <= tsplice))
  # t^l < t < x_i=l_i=u_i <=T
  ind2 <- which(uncensored & (upper > tsplice))
  # t^l <= l_i < u_i <= t <T
  ind3 <- which(censored & upper <= tsplice)
  # t^l < t <= l_i < u_i <=T
  ind4 <- which(censored & lower >= tsplice)
  # t^l <= l_i < t < u_i <=T
  ind5 <- which(censored & lower < tsplice & upper > tsplice)

  lower1 <- lower[ind1];
  lower2 <- lower[ind2];
  lower3 <- lower[ind3]; upper3 <- upper[ind3]
  lower4 <- lower[ind4]; upper4 <- upper[ind4]
  lower5 <- lower[ind5]; upper5 <- upper[ind5]

  # Call C++ code
  .Call("_ReIns_spliceEM_splicefit_raw_Rexport", pi, theta, shape, beta, gamma,
        lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5, 
        trunclower, tsplice, truncupper, eps, beta_tol, maxiter,
        PACKAGE = "ReIns")
}




.spliceEM_shape_adj_cpp <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf,
                                    pi, theta, shape, beta, gamma,
                                    eps = 10 ^ (-3), beta_tol = 10 ^ (-5), maxiter = Inf) {

  # Separate uncensored and censored observations
  uncensored <- !censored

  # Indices for the five different cases
  # Some can be numeric(0)!

  # t^l <= x_i=l_i=u_i <= t <T
  ind1 <- which(uncensored & (upper <= tsplice))
  # t^l < t < x_i=l_i=u_i <=T
  ind2 <- which(uncensored & (upper > tsplice))
  # t^l <= l_i < u_i <= t <T
  ind3 <- which(censored & upper <= tsplice)
  # t^l < t <= l_i < u_i <=T
  ind4 <- which(censored & lower >= tsplice)
  # t^l <= l_i < t < u_i <=T
  ind5 <- which(censored & lower < tsplice & upper > tsplice)

  lower1 <- lower[ind1];
  lower2 <- lower[ind2];
  lower3 <- lower[ind3]; upper3 <- upper[ind3]
  lower4 <- lower[ind4]; upper4 <- upper[ind4]
  lower5 <- lower[ind5]; upper5 <- upper[ind5]

  # Call C++ code
  .Call("_ReIns_spliceEM_shape_adj_Rexport", pi, theta, shape, beta, gamma, 
        lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
        trunclower, tsplice, truncupper, eps, beta_tol, maxiter,
        PACKAGE = "ReIns")
}


.spliceEM_shape_red_cpp <- function(lower, upper, censored, trunclower = 0, tsplice, truncupper = Inf, 
                                   pi, theta, shape, beta, gamma, criterium = "AIC", improve = TRUE,
                                   eps = 10 ^ (-3), beta_tol = 10 ^ (-5), maxiter = Inf, adj = TRUE) {
  
  # Separate uncensored and censored observations
  uncensored <- !censored
  
  # Indices for the five different cases
  # Some can be numeric(0)!
  
  # t^l <= x_i=l_i=u_i <= t <T
  ind1 <- which(uncensored & (upper <= tsplice))
  # t^l < t < x_i=l_i=u_i <=T
  ind2 <- which(uncensored & (upper > tsplice))
  # t^l <= l_i < u_i <= t <T
  ind3 <- which(censored & upper <= tsplice)
  # t^l < t <= l_i < u_i <=T
  ind4 <- which(censored & lower >= tsplice)
  # t^l <= l_i < t < u_i <=T
  ind5 <- which(censored & lower < tsplice & upper > tsplice)
  
  lower1 <- lower[ind1];
  lower2 <- lower[ind2]; 
  lower3 <- lower[ind3]; upper3 <- upper[ind3]
  lower4 <- lower[ind4]; upper4 <- upper[ind4]
  lower5 <- lower[ind5]; upper5 <- upper[ind5]
  
  # Call C++ code
  .Call("_ReIns_spliceEM_shape_red", pi, theta, shape, beta, gamma, 
        lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5, 
        trunclower, tsplice, truncupper, criterium, improve, eps, beta_tol, maxiter, adj,
        PACKAGE = "ReIns")
}


# Unload package DLL when package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("ReIns", libpath)
}