
#include "SpliceEM_Estep.h"

// Compute P(X_i<=t | t^l<=l_i<=t<u_i, Theta^{(h-1)}) and P(X_i>t | t^l<=l_i<=t<u_i, Theta^{(h-1)})
// Use pi^{(h-1)}, gamma^{(h-1)}, ...
NumericVector spliceEM_probs(const NumericVector &lower5, const NumericVector &upper5, 
                             const double trunclower, const double tsplice, const double truncupper, 
                             const double pi, const double theta, const NumericVector shape, const NumericVector alpha, const double gamma) {

  // CDF of ME evaluated in l_i (for all i with t^l<=l_i<=t<u_i)
  NumericVector p1 = pME_vec(lower5, theta, shape, alpha, trunclower, tsplice);
  
  NumericVector P1 = pi*(1-p1) / (pi+(1-pi)*ptpareto_vec(upper5, gamma, tsplice, truncupper) - pi*p1);

  return P1;
}



// ^{i.}z_{ij}^{(h)}: posterior probabilities (uncensored)
NumericMatrix spliceEM_i_z(const NumericMatrix &x1_dens_nosum, const int M) {

  NumericMatrix z1 = x1_dens_nosum; 

  for (int i = 0; i < z1.nrow(); ++i) {
    
    z1.row(i) = z1.row(i) / sum(z1.row(i));
    
    // in case all ^{u}z_{ij}^{(h)} for j=1,...,M are numerically 0
    if (sum(Rcpp::is_nan(z1.row(i))) == M) {
      z1.row(i) = rep(1/M, M);
    }
    
  } 

  return z1;
}

// ^{iii.}z_{ij}^{(h)}: posterior probabilities (censored)
NumericMatrix spliceEM_iii_z(const NumericMatrix &c3_probs_nosum, const int M) {

  NumericMatrix z3 = c3_probs_nosum; 
  
  for (int i = 0; i < z3.nrow(); ++i) {
    
    z3.row(i) = z3.row(i) / sum(z3.row(i));
    
    // in case all ^{c}z_{ij}^{(h)} for j=1,...,M are numerically 0
    if (sum(Rcpp::is_nan(z3.row(i))) == M) {
      z3.row(i) = rep(1/M, M);
    }
  } 
  return z3;
}


// ^{v.}z_{ij}^{(h)}: posterior probabilities (censored)
NumericMatrix spliceEM_v_z(const NumericMatrix &c5_probs_nosum, const double tsplice, const NumericVector alpha_tilde, 
                           const NumericVector shape, const double theta,  const int M) {

  NumericVector pt(M);
  
  for (int j = 0; j < M; ++j) {
    // Use alpha_tilde since also used in c5_probs_nosum
    pt[j] = alpha_tilde[j] * pGamma(tsplice, shape[j], theta);
    
  }

  NumericMatrix z5 = c5_probs_nosum;
  
  for (int i = 0; i < z5.nrow(); ++i) {
    
    z5.row(i) = pt - z5.row(i);
    
    z5.row(i) = z5.row(i) / sum(z5.row(i));
    
    // in case all ^{cc}z_{ij}^{(h)} for j=1,...,M are numerically 0
    if (sum(Rcpp::is_nan(z5.row(i))) == M) {
      z5.row(i) = rep(1/M, M);
    }
    
  } 

  return z5;
}


// Expected value of censored observations with l_i<u_i<t for ME
NumericMatrix spliceEM_Estep_ME_iii(const NumericVector &lower3, const NumericVector &upper3, 
                                    const NumericVector shape, const double theta) {

  int n3 = lower3.size();
  int M = shape.size();
  NumericMatrix E3_ME = NumericMatrix(n3,M);
  
  for (int i = 0; i < n3; ++i) {
    
    for (int j = 0; j < M; ++j) {
  
      E3_ME(i,j) = shape[j] * theta * (pGamma(upper3[i], shape[j]+1, theta) - pGamma(lower3[i], shape[j]+1, theta)) / (pGamma(upper3[i], shape[j], theta) - pGamma(lower3[i], shape[j], theta));
    
      // replace numerical 0/0 (NaN) or Inf by correct expected value
      if (::R_IsNaN(E3_ME(i,j)) || !::R_finite(E3_ME(i,j))) {
        
        if (lower3[i] > shape[j]*theta) {
          E3_ME(i,j) = lower3[i];
          
        } else {
          E3_ME(i,j) = upper3[i];
        }
      }
      
    }
  }

  return E3_ME;
}

// Expected value of censored observations with l_i<t<u_i for ME
NumericMatrix spliceEM_Estep_ME_v(const NumericVector &lower5, const double tsplice, 
                                  const NumericVector shape, const double theta) {

  return spliceEM_Estep_ME_iii(lower5, rep(tsplice, lower5.size()), shape, theta);

}


// E-step for Pareto: part of E(ln f_2(X_i;t) | t^l<t<=l_i<u_i, Theta^{(h-1)}) not depending on gamma
NumericVector spliceEM_Estep_Pa_iv(const NumericVector &lower4, const NumericVector &upper4, 
                                   const double gamma, const double tsplice) {

  // Problem if upper4 is infinite
  NumericVector u4(upper4.size());
  NumericVector u4_gamma(upper4.size());
  
  for (int i = 0; i < upper4.size(); ++i) {

    if (::R_finite(upper4[i])) {
      
      u4[i] = upper4[i]/tsplice;
      u4_gamma[i] = pow(upper4[i]/tsplice, -1/gamma);
      
    } else {
      
      u4[i] = 1;
      u4_gamma[i] = 0;
    }
  }
  
  return ((log(lower4/tsplice)+gamma) * pow(lower4/tsplice, -1/gamma) - (log(u4)+gamma) * u4_gamma) / (pow(lower4/tsplice, -1/gamma) - u4_gamma);
}


// E-step for Pareto: part of E(ln f_2(X_i;t) | t^l<l_i<=t<u_i, Theta^{(h-1)}) not depending on gamma
NumericVector spliceEM_Estep_Pa_v(const NumericVector upper5, const double gamma, const double tsplice) {

  return spliceEM_Estep_Pa_iv(rep(tsplice, upper5.size()), upper5, gamma, tsplice);

}