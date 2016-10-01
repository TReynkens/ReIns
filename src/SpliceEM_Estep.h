#ifndef SpliceEM_Estep_H
#define SpliceEM_Estep_H

#include "SpliceEM_aux.h"

// Compute P(X_i<=t | t^l<=l_i<=t<u_i, Theta^{(h-1)}) and P(X_i>t | t^l<=l_i<=t<u_i, Theta^{(h-1)})
NumericVector spliceEM_probs(const NumericVector &lower5, const NumericVector &upper5, 
                             const double trunclower, const double tsplice, const double truncupper, 
                             const double pi, const double theta, const NumericVector shape, const NumericVector alpha, const double gamma);


// ^{i.}z_{ij}^{(h)}: posterior probabilities (uncensored)
NumericMatrix spliceEM_i_z(const NumericMatrix &x1_dens_nosum, const int M);
  
// ^{iii.}z_{ij}^{(h)}: posterior probabilities (censored)  
NumericMatrix spliceEM_iii_z(const NumericMatrix &c3_probs_nosum, const int M);
  
// ^{v.}z_{ij}^{(h)}: posterior probabilities (censored)  
NumericMatrix spliceEM_v_z(const NumericMatrix &c5_probs_nosum, const double tsplice, const NumericVector alpha_tilde, 
                           const NumericVector shape, const double theta, const int M);

// Expected value of censored observations with l_i<u_i<t for ME
NumericMatrix spliceEM_Estep_ME_iii(const NumericVector &lower3, const NumericVector &upper3, 
                                    const NumericVector shape, const double theta);
  
// Expected value of censored observations with l_i<t<u_i for ME  
NumericMatrix spliceEM_Estep_ME_v(const NumericVector &lower5, const double tsplice, 
                                  const NumericVector shape, const double theta);


// E-step for Pareto: part of E(ln f_2(X_i;t) | t^l<t<=l_i<u_i, Theta^{(h-1)}) not depending on gamma  
NumericVector spliceEM_Estep_Pa_iv(const NumericVector &lower4, const NumericVector &upper4, 
                                   const double gamma, const double tsplice);

// E-step for Pareto: part of E(ln f_2(X_i;t) | t^l<l_i<=t<u_i, Theta^{(h-1)}) not depending on gamma    
NumericVector spliceEM_Estep_Pa_v(const NumericVector &upper5, const double gamma, const double tsplice);
      
#endif