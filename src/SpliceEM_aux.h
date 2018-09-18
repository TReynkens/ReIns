#ifndef SpliceEM_aux_H // Check if h-file is defined before
#define SpliceEM_aux_H // Define h-file

// Include libraries and header files
#include <iostream> // Input and output
#include <cmath> // C numerics library
#include <cfloat> // Characteristics of floating-point types

// R functionality
#include <Rcpp.h> // R/C++ interface class library
using namespace Rcpp;


///////////////////////////////////////////////////////
// Matrix calculations

// Row and column sums of matix
NumericVector colSums(const NumericMatrix &X);
NumericVector rowSums(const NumericMatrix &X);


//////////////////////////////////////////
// Distributions

// Gamma distribution (PDF and CDF)
double dGamma(const double x, const double shape, const double scale);
double pGamma(const double x, const double shape, const double scale);


// PDF of Pareto distribution and vector input version
double dpareto(const double x, const double gamma, const double t);
NumericVector dpareto_vec(const NumericVector x, const double gamma, const double t);

// CDF of Pareto distribution and vector input version
double ppareto(const double x, const double gamma, const double t);  
NumericVector ppareto_vec(const NumericVector x, const double gamma, const double t);


// PDF of truncated Pareto distribution and vector input version
double dtpareto(const double x, const double gamma, const double t, const double truncupper);
NumericVector dtpareto_vec(const NumericVector x, const double gamma, const double t, const double truncupper);

// CDF of truncated Pareto distribution and vector input version
double ptpareto(const double x, const double gamma, const double t, const double truncupper);
NumericVector ptpareto_vec(const NumericVector x, const double gamma, const double t, const double truncupper);


// CDF of Mixed Erlang distribution and vector input version
double pME(const double x, const double theta, const IntegerVector shape, 
           const NumericVector alpha, const double trunclower, const double truncupper);
NumericVector pME_vec(const NumericVector x, const double theta, const IntegerVector shape, 
                      const NumericVector alpha, const double trunclower, const double truncupper);

#endif