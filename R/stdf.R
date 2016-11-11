########################
# Estimators of the stable tail dependence function

# Call C++ and make sure the right input types are used
stdf <- function(x, k, X, alpha=0.5) {

  # Convert to rank matrix
  R <- apply(as.matrix(X), 2, rank)
  return(.Call("ReIns_stdf", as.numeric(x), k, as.matrix(R), alpha))
}

# Call C++ and make sure the right input types are used
stdf2 <- function(x, k, X) {

  return(.Call("ReIns_stdf2", as.numeric(x), k, as.matrix(X)))
}


