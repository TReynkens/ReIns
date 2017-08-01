########################
# Estimators of the stable tail dependence function

# Call C++ and make sure the right input types are used
stdf <- function(x, k, X, alpha = 0.5) {

  # Convert to rank matrix
  R <- apply(as.matrix(X), 2, rank)
  return(.Call("_ReIns_stdf_cpp", as.numeric(x), k, as.matrix(R), alpha,
               PACKAGE = "ReIns"))
}

# Call C++ and make sure the right input types are used
stdf2 <- function(x, k, X) {

  return(.Call("_ReIns_stdf2_cpp", as.numeric(x), k, as.matrix(X),
               PACKAGE = "ReIns"))
}
