

// Include libraries and header files
#include <iostream> // Input and output
#include <cmath> // C numerics library

// R functionality
#include <Rcpp.h> // R/C++ interface class library
using namespace Rcpp;


///////////////////////////////////////////////////////

// Estimator for the stable tail dependence function: hat(l)
// R is the matrix of ranks!
// [[Rcpp::export]]
double stdf_cpp(const NumericVector &x, const int k, const NumericMatrix &R, const double alpha) {

  int n = R.nrow();
  int d = R.ncol();
  
  // Check if x has the right length
  if (x.size()!=d) {
    ::Rf_error("x should be a vector with the same length as the number of columns of R.");
  } 
  
  // Compute thresholds to use
  NumericVector threshold = n + alpha - k * x;
  
  double l = 0;
  // Look in each row if at least one element above (or equal to) threshold
  for (int i = 0; i < n; ++i) {
    l += sum(R.row(i)<=threshold)!=d;
  }
  return l / k;
}


// Estimator for the stable tail dependence function: tilde(l)
// [[Rcpp::export]]
double stdf2_cpp(const NumericVector &x, const int k, NumericMatrix &X) {
  
  int n = X.nrow();
  int d = X.ncol();
  
  // Check if x has the right length
  if (x.size()!=d) {
    ::Rf_error("x should be a vector with the same length as the number of columns of X.");
  } 
  
  // Compute thresholds to use
  NumericVector threshold(d);
  
  for (int j = 0; j < d; ++j) {
    // Select j-th row and do not copy it
    NumericVector Xrow = X(_,j);
    // Index of element to look at
    int ind = (n-ceil(k*x[j])+1)-1;
    // Avoid negative indices
    if (ind<0) ind = 0;
    // Partially sort vector
    std::nth_element(Xrow.begin(), Xrow.begin() + ind, Xrow.end());
    
    // Select right element as threshold
    threshold[j] = Xrow[ind];
  }
  
  
  double l = 0;
  // Look in each row if at least one element above (or equal to) threshold
  for (int i = 0; i < n; ++i) {
    l += sum(X.row(i)<threshold)!=d;
  }
  return l / k;
}

