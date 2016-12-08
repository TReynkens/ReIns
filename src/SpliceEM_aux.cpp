// Include header file
#include "SpliceEM_aux.h"


///////////////////////////////////////////////////////
// Matrix calculations

// Column sums of a matrix
NumericVector colSums(const NumericMatrix &X) {
  
  int ncol = X.ncol();
  
  NumericVector colsums(ncol); 
  
  for (int j = 0; j < ncol; ++j) {
    colsums[j] = sum(X.column(j));
  }
  
  return colsums;
}

// Row sums of a matrix
NumericVector rowSums(const NumericMatrix &X) {
  
  int nrow = X.nrow();
  
  NumericVector rowsums(nrow); 
  
  for (int i = 0; i < nrow; ++i) {
    rowsums[i] = sum(X.row(i));
  }
  
  return rowsums;
}


//////////////////////////////////////////
// Distributions


// Gamma PDF
double dGamma(const double x, const double shape, const double scale) {
  
  // Gamma PDF with arguments shape and scale
  // Use 0 for log.p
  return R::dgamma(x, shape, scale, 0);
}


// Gamma CDF
double pGamma(const double x, const double shape, const double scale) {
  
  // Gamma CDF with arguments shape and scale
  // Use 1 for lower.tail and 0 for log.p
  return R::pgamma(x, shape, scale, 1, 0);
}



// Pareto PDF
double dpareto(const double x, const double gamma, const double t) {
  
  double d;
  
  if (x <= t) {
    d = 0;
    
  } else  {
    d = 1 / (gamma*t) * pow(x/t, -1/gamma-1);
  }
  
  return d;
}

// Pareto PDF (vector input)
NumericVector dpareto_vec(const NumericVector x, const double gamma, const double t) {
  
  int n = x.size();
  
  NumericVector d(n);
  
  for (int i = 0; i < n; ++i) {
    
    d[i] =  dpareto(x[i], gamma, t);
    
  }
  
  return d;
} 


// Pareto CDF
double ppareto(const double x, const double gamma, const double t) {
  
  double p;
  
  if (x <= t) {
    p = 0;
    
  } else  {
    p = 1 - pow(x/t, -1/gamma);
  }
  
  return p;
} 

// Pareto CDF (vector version)
NumericVector ppareto_vec(const NumericVector x, const double gamma, const double t) {
  
  int n = x.size();
  
  NumericVector p(n);
  
  for (int i = 0; i < n; ++i) {
    
    p[i] = ppareto(x[i], gamma, t);
    
  }
  
  return p;
} 



// Truncated Pareto PDF
double dtpareto(const double x, const double gamma, const double t, const double truncupper) {
  
  double d;
  
  if (x <= t || x >= truncupper) {
    d = 0;
    
  } else  {
    d = dpareto(x, gamma, t) / ppareto(truncupper, gamma, t);
  }
  
  return d;
} 

// Truncated Pareto PDF (vector version)
NumericVector dtpareto_vec(const NumericVector x, const double gamma, const double t, const double truncupper) {
  
  int n = x.size();
  
  NumericVector d(n);
  
  for (int i = 0; i < n; ++i) {
    
    d[i] = dtpareto(x[i], gamma, t, truncupper);
  }
  
  return d;
} 


// Truncated Pareto CDF
double ptpareto(const double x, const double gamma, const double t, const double truncupper) {
  
  double p;
  
  if (x <= t) {
    p = 0;
    
  } else if(x < truncupper) {
    p = ppareto(x, gamma, t) / ppareto(truncupper, gamma, t);
    
  } else {
    p = 1;
  }
  
  return p;
} 

// Truncated Pareto CDF (vector version)
NumericVector ptpareto_vec(const NumericVector x, const double gamma, const double t, const double truncupper) {
  
  int n = x.size();
  
  NumericVector p(n);
  
  for (int i = 0; i < n; ++i) {
    
    p[i] = ptpareto(x[i], gamma, t, truncupper);
  }
  
  return p;
} 



// Mixed Erlang CDF
double pME(const double x, const double theta, const IntegerVector shape, 
           const NumericVector alpha, const double trunclower, const double truncupper) {
  
  double p = 0;
  
  int M = shape.size();
  
  // Compute CDF as weighted sum of gamma CDF's
  for(int j = 0; j < M; j++) {
    
    p += pGamma(x, shape[j], theta) * alpha[j];
  }
  
  
  // Truncation
  if (trunclower!= 0 || R_FINITE(truncupper)) {
    
    double pl = pME(trunclower, theta, shape, alpha, 0, R_PosInf);
    double pu = pME(truncupper, theta, shape, alpha, 0, R_PosInf);
    
    if (x <= trunclower) {
      p = 0;
      
    } else if (x < truncupper) {
      p = (p - pl) / (pu - pl);
      
    } else {
      p = 1;
    }
  }
  
  return p;
}

// Mixed Erlang CDF (vector version)
NumericVector pME_vec(const NumericVector x, const double theta, const IntegerVector shape, 
                      const NumericVector alpha, const double trunclower, const double truncupper) {
  
  int n = x.size();
  NumericVector p(n);
  
  for(int i = 0; i < n; i++) {
    
    p[i] = pME(x[i], theta, shape, alpha, trunclower, truncupper);
  }
  
  return p;
}

