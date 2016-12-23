
#include "SpliceEM_aux.h"
#include "SpliceEM_Estep.h"
#include "SpliceEM_Mstep.h"

// Compute splicing log-likelihood
double splice_loglikelihood(const NumericVector &x1_dens, const NumericVector &x2_dens, const NumericVector &c3_probs, 
                            const NumericVector &c4_probs, const NumericVector &c5_probs) {
  
  int n1 = x1_dens.size();
  int n2 = x2_dens.size();
  int n3 = c3_probs.size();
  int n4 = c4_probs.size();
  int n5 = c5_probs.size();
  
  double ll = 0;
  
  
  // Determine log-likelihood contributions for all five classes
  // The log of a negative contribution is replaced by -1000
  
  // Contributions of observations of classes 1 and 2 are log-densities
  if (n1>0) {
    
    for (int i = 0; i < n1; ++i) {
      
      if (x1_dens[i]>0) {
        
        ll += log(x1_dens[i]);
        
      } else {
        
        ll += -1000;
      }
    }
  }
  
  
  if (n2>0) {
    
    for (int i = 0; i < n2; ++i) {
      
      if (x2_dens[i]>0) {
        
        ll += log(x2_dens[i]);
        
      } else {
        
        ll += -1000;
      }
    }
  }
  
  
  // Contributions of observations of classes 3, 4 and 5 are log-probabilities (censoring)
  if (n3>0) {
    
    for (int i = 0; i < n3; ++i) {
      
      if (c3_probs[i]>0) {
        
        ll += log(c3_probs[i]);
        
      } else {
        
        ll += -1000;
      }
    }
  }
  
  if (n4>0) {
    
    for (int i = 0; i < n4; ++i) {
      
      if (c4_probs[i]>0) {
        
        ll += log(c4_probs[i]);
        
      } else {
        
        ll += -1000;
      }
    }
  }
  
  
  if (n5>0) {
    
    for (int i = 0; i < n5; ++i) {
      
      if (c5_probs[i]>0) {
        
        ll += log(c5_probs[i]);
        
      } else {
        
        ll += -1000;
      }
    }
  }
  
  
  return ll; 
}


// Compute densities and probabilities for observations from the 5 classes
void spliceEM_densprob(NumericMatrix &x1_dens_nosum, NumericVector &x1_dens, NumericVector &x2_dens, NumericMatrix &c3_probs_nosum, NumericVector &c3_probs,
                       NumericVector &c4_probs, NumericMatrix &c5_probs_nosum, NumericVector &c5_probs,
                       const NumericVector &lower1, const NumericVector &lower2, const NumericVector &lower3, const NumericVector &lower4, const NumericVector &lower5,
                       const NumericVector &upper3, const NumericVector &upper4, const NumericVector &upper5,
                       const double trunclower, const double tsplice, const double truncupper,
                       const double pi, const double theta, const IntegerVector shape, const NumericVector alpha_tilde, const double gamma) {
  
  int M = shape.size();
  int n1 = lower1.size();
  int n3 = lower3.size();
  int n5 = lower5.size();
  
  if (n1>0) {
    
    // Matrix containing Erlang densities for uncensored observations in ME part and all M values of shape
    // Note F_1(x) = sum_{j=1}^M beta_j F_E^t(x) = sum_{j=1}^M alpha_j F_E(x)/(F_E(t)-F_E(t^l))
    // and alpha_tilde = beta/(F_E(t)-F_E(t^l))
    x1_dens_nosum = NumericMatrix(n1, M);
    
    for (int i = 0; i < n1; ++i) {
      
      for (int j = 0; j < M; ++j) {
        x1_dens_nosum(i,j) = dGamma(lower1[i], shape[j], theta) * alpha_tilde[j];
      }
    }
    
    // Combine into density vector (multiply with pi because splicing density)
    x1_dens = pi * rowSums(x1_dens_nosum);
  }
  
  
  // Densities of uncensored observations in EVT part
  x2_dens = (1-pi) * dtpareto_vec(lower2, gamma, tsplice, truncupper);
  
  
  if (n3>0) {
    
    // Matrix containing Erlang probabilities for censored observations of type 3 and all M values of shape
    c3_probs_nosum = NumericMatrix(n3, M);
    
    for (int i = 0; i < n3; ++i) {
      
      for (int j = 0; j < M; ++j) {
        c3_probs_nosum(i,j) = (pGamma(upper3[i], shape[j], theta) - pGamma(lower3[i], shape[j], theta)) * alpha_tilde[j];
      }
      
    }
    // Combine into density vector (multiply with pi because splicing density)
    c3_probs = pi * rowSums(c3_probs_nosum);
  }
  
  
  // Vector of probabilities for uncensored observations of type 4: F(u)-F(l)
  c4_probs = pi + (1-pi) * (ptpareto_vec(upper4, gamma, tsplice, truncupper) - ptpareto_vec(lower4, gamma, tsplice, truncupper));
  
  
  if (n5>0) {
    
    // Matrix containing Erlang probabilities for censored observations of type 5 and all M values of shape
    c5_probs_nosum = NumericMatrix(n5, M);
    
    for (int i = 0; i < n5; ++i) {
      
      for (int j = 0; j < M; ++j) {
        c5_probs_nosum(i,j) = pGamma(lower5[i], shape[j], theta) * alpha_tilde[j];
      }
      
    }
    // Vector of probabilities for uncensored observations of type 5: F(u)-F(l)
    c5_probs = pi + (1-pi) * ptpareto_vec(upper5, gamma, tsplice, truncupper) - pi * rowSums(c5_probs_nosum);
  }
  
}




// Fit splicing model using EM algorithm for given M and shapes
// Additional input: eps = numerical tolerance for EM algorithm (stopping criterion),
// beta_tol = numerical tolerence to discard beta's,
// maxiter = maximal number of EM iterations
List spliceEM_splicefit_raw(const double pi_in, const double theta_in, const IntegerVector shape_in, const NumericVector beta_in, const double gamma_in,
                            const NumericVector &lower1, const NumericVector &lower2, const NumericVector &lower3, const NumericVector &lower4, const NumericVector &lower5, 
                            const NumericVector &upper3, const NumericVector &upper4, const NumericVector &upper5, 
                            const double trunclower, const double tsplice, const double truncupper,
                            const double eps, const double beta_tol, const double maxiter) {
  
  double pi = pi_in;
  double theta = theta_in;
  // Clone vectors to avoid changing it in R (or C++)
  IntegerVector shape = clone(shape_in);
  NumericVector beta = clone(beta_in);
  double gamma = gamma_in;
  
  int M = beta.size();
  // n1 and n2 have different meaning in paper!
  int n1 = lower1.size();
  int n2 = lower2.size();
  int n3 = lower3.size();
  int n4 = lower4.size();
  int n5 = lower5.size();
  int n = n1 + n2 + n3 + n4 + n5;
  
  // Truncation probabilities
  NumericVector t_probs(M);
  
  for (int j = 0; j < M; ++j) {
    
    t_probs[j] = pGamma(tsplice, shape[j], theta) - pGamma(trunclower, shape[j], theta);
  }
  
  
  
  // MLE for alpha, transform beta to alpha
  // alpha_tilde = alpha/(F_ME(t)-F_ME(t^l))
  NumericVector alpha_tilde = beta / t_probs;
  NumericVector alpha = alpha_tilde / sum(alpha_tilde);
  
  // Compute densities and probabilities for observations from the 5 classes
  NumericMatrix x1_dens_nosum;
  NumericVector x1_dens;
  
  NumericVector x2_dens;
  
  NumericMatrix c3_probs_nosum;
  NumericVector c3_probs;
  
  NumericVector c4_probs;
  
  NumericMatrix c5_probs_nosum;
  NumericVector c5_probs;
  
  
  // Compute densities and probabilities for observations from the 5 classes
  spliceEM_densprob(x1_dens_nosum, x1_dens, x2_dens, c3_probs_nosum, c3_probs, c4_probs, c5_probs_nosum, c5_probs,
                    lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                    trunclower, tsplice, truncupper, pi, theta, shape, alpha_tilde, gamma);
  

  // Compute log-likelihood
  double loglikelihood = splice_loglikelihood(x1_dens, x2_dens, c3_probs, c4_probs, c5_probs);
  
  int iteration = 1;
  double old_loglikelihood = R_NegInf;
  
  // Start EM-loop
  
  while((loglikelihood - old_loglikelihood > eps) && (iteration <= maxiter)) {
    
    old_loglikelihood = loglikelihood;
    
    ////////////////
    // E-step
    ////////////////
    
    NumericVector P1, P2;
    
    if (n5>0) {
      
      // Compute P(X_i<=t | t^l<=l_i<=t<u_i, Theta^{(h-1)}) and P(X_i>t | t^l<=l_i<=t<u_i, Theta^{(h-1)})
      P1 = spliceEM_probs(lower5, upper5, trunclower, tsplice, truncupper, pi, theta, shape, alpha, gamma);
      
      P2 = 1 - P1;
      
    } else {
      
      P1 = 0, P2 = 0;
      
    }
    
    ///////
    // pi
    
    double n1_h = n1 + n3 + sum(P1);
    double n2_h = n - n1_h;
    
    ///////
    // ME
    
    NumericMatrix z1(1, M), z3(1, M), z5(1, M);
    double E1_ME, E3_ME, E5_ME;
    
    E1_ME = 0;
    if (n1>0) {
      z1 = spliceEM_i_z(x1_dens_nosum, M);
      
      // sum_{i in S_{i.}} x_i
      E1_ME = sum(lower1);
      
    } 
    
    E3_ME = 0;
    if (n3>0) {
      z3 = spliceEM_iii_z(c3_probs_nosum, M);
      
      // E(X_{ij} | Z_{ij}=1, t^l<=l_i<u_i<=t)
      NumericMatrix Estep_ME_iii = spliceEM_Estep_ME_iii(lower3, upper3, shape, theta);
      
      // sum_{i in S_{iii.}} sum_{j=1}^m z3 * E(X_{ij} | Z_{ij}=1, t^l<=l_i<u_i<=t)
      for (int i = 0; i < n3; ++i) {
        
        E3_ME += sum(z3.row(i) * Estep_ME_iii.row(i));
      }
      
    } 

    
    E5_ME = 0;
    if (n5>0) {
      z5 = spliceEM_v_z(c5_probs_nosum, tsplice, alpha_tilde, shape, theta, M);
      
      // E(X_{ij} | Z_{ij}=1, t^l<=l_i<t<u_i)
      NumericMatrix Estep_ME_v = spliceEM_Estep_ME_v(lower5, tsplice, shape, theta);
      
      // sum_{i in S_{v.}} P_{1,i} * sum_{j=1}^m z3 * E(X_{ij} | Z_{ij}=1, t^l<=l_i<t<u_i)
      for (int i = 0; i < n5; ++i) {
        
        E5_ME += P1[i] * sum(z5.row(i) * Estep_ME_v.row(i));
      }
    }


    ///////
    // gamma
    
    double E2_Pa, E4_Pa, E5_Pa;
    
    if (n2>0) {
      // sum_{i in S_{ii.}} ln(x_i/t)
      E2_Pa = sum(log(lower2/tsplice));
      
    } else {
      E2_Pa = 0;
    }
    
    if (n4>0) {
      // sum_{i in S_{iv.}} E( ln(X_i/t) | t^l<t<=l_i<u_i, gamma^{(h-1)})
      E4_Pa = sum(spliceEM_Estep_Pa_iv(lower4, upper4, gamma, tsplice));
      
    } else {
      E4_Pa = 0;
    }
    
    if (n5>0) {
      // sum_{i in S_{v.}} P_{2,i} * E( ln(X_i/t) | t^l<l_i<=t<u_i, gamma^{(h-1)})
      E5_Pa = sum(P2 * spliceEM_Estep_Pa_v(upper5, gamma, tsplice));
      
    } else {
      E5_Pa = 0;
    }
    
    

    ////////////////
    // M-step
    ////////////////
    
    
    
    ///////
    // pi
    
    pi = n1_h / n;
    
    ///////
    // ME

    for (int j = 0; j < M; ++j) {
      
      beta[j] = sum(z1.column(j)) + sum(z3.column(j)) + sum(z5.column(j)*P1);
    }
    
    beta = beta / n1_h;
    
    // Remove small beta's
    if (min(beta) < beta_tol) {
      shape = shape[beta > beta_tol];
      beta = beta[beta > beta_tol];
      beta = beta / sum(beta);
      // New M
      M = beta.size();
    }
    
    
    theta = exp(spliceEM_theta(log(theta), E1_ME, E3_ME, E5_ME, n1_h, beta, shape, trunclower, tsplice));

    ///////
    // gamma  
    
    double gamma_old = gamma;
    
    gamma = 1/n2_h * (E2_Pa + E4_Pa + E5_Pa);
    
    
    // Solve equation numerically 
    if (R_FINITE(truncupper)) {
      // spliceEM_Mstep_Pareto gives back log(gamma_opt)
      gamma = exp(spliceEM_Mstep_Pareto(log(gamma_old), gamma, tsplice, truncupper));
    } 

    ////////////////
    
    ++iteration;
    
    // Truncation probabilities
    for (int j = 0; j < M; ++j) {
      
      t_probs[j] = pGamma(tsplice, shape[j], theta) - pGamma(trunclower, shape[j], theta);
    }
    
    // MLE for alpha, transform beta to alpha
    // alpha_tilde = alpha/(F_ME(t)-F_ME(t^l))
    alpha_tilde = beta / t_probs;
    alpha = alpha_tilde / sum(alpha_tilde);
    
    // Compute densities and probabilities for observations from the 5 classes
    spliceEM_densprob(x1_dens_nosum, x1_dens, x2_dens, c3_probs_nosum, c3_probs, c4_probs, c5_probs_nosum, c5_probs,
                      lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                      trunclower, tsplice, truncupper, pi, theta, shape, alpha_tilde, gamma);

    // Compute log-likelihood
    loglikelihood = splice_loglikelihood(x1_dens, x2_dens, c3_probs, c4_probs, c5_probs);
    
  }
  
  // Compute ICs
  int df = 2*alpha.size()-1+1+1+1;
  double aic = -2 * loglikelihood + 2.0 * df;
  double bic = -2 * loglikelihood + std::log(n*1.0) * df;
  
  
  // Make output list
  List output = Rcpp::List::create(Rcpp::Named("pi") = pi,
                                   Rcpp::Named("alpha") = alpha,
                                   Rcpp::Named("beta") = beta,
                                   Rcpp::Named("shape") = shape,
                                   Rcpp::Named("theta") = theta,
                                   Rcpp::Named("gamma") = gamma,
                                   Rcpp::Named("loglikelihood") = loglikelihood,
                                   Rcpp::Named("iteration") = iteration-1,
                                   Rcpp::Named("AIC") = aic,
                                   Rcpp::Named("BIC") = bic);

  return output;
}

// Fit splicing model using EM algorithm for given M and shapes (export to R)
// [[Rcpp::export]]
List spliceEM_splicefit_raw_Rexport(const double pi, const double theta, const IntegerVector shape, const NumericVector beta, const double gamma, 
                                    const NumericVector lower1, const NumericVector lower2, const NumericVector lower3, const NumericVector lower4, const NumericVector lower5, 
                                    const NumericVector upper3, const NumericVector upper4, const NumericVector upper5, 
                                    const double trunclower, const double tsplice, const double truncupper,
                                    const double eps, const double beta_tol, const double maxiter) {
  
  
  // Call fitting function
  return spliceEM_splicefit_raw(pi, theta, shape, beta, gamma,
                                lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
  
}



// Shape adjustments (4.2 in Verbelen et al. 2015)
List spliceEM_shape_adj(const double pi_in, const double theta_in, const IntegerVector shape_in, const NumericVector beta_in, const double gamma_in, 
                        const NumericVector &lower1, const NumericVector &lower2, const NumericVector &lower3, const NumericVector &lower4, const NumericVector &lower5, 
                        const NumericVector &upper3, const NumericVector &upper4, const NumericVector &upper5, 
                        const double trunclower, const double tsplice, const double truncupper,
                        const double eps, const double beta_tol, const double maxiter) {
  
  double pi = pi_in;
  double theta = theta_in;
  // Clone vectors to avoid changing it in R (and C++)
  IntegerVector shape = clone(shape_in);
  NumericVector beta = clone(beta_in);
  double gamma = gamma_in;
  
  // Fit model with shape equal to the initial value (and M=length(shape))
  List fit = spliceEM_splicefit_raw(pi, theta, shape, beta, gamma,
                                    lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                    trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
  
  double loglikelihood = as<double>(fit["loglikelihood"]);
  pi = as<double>(fit["pi"]);
  theta = as<double>(fit["theta"]);
  shape = as<IntegerVector>(fit["shape"]);
  beta = as<NumericVector>(fit["beta"]);
  NumericVector alpha = as<NumericVector>(fit["alpha"]);
  gamma = as<double>(fit["gamma"]);
  int M = beta.size();

  int iteration = 1;
  double before_loglikelihood = R_NegInf;
  double after_loglikelihood = loglikelihood;
  double new_loglikelihood;
  
  // Improve as long as log-likelihood increases significantly
  while(after_loglikelihood - before_loglikelihood > eps) {
    
    before_loglikelihood = after_loglikelihood;
    
    // Try increasing the shapes (step 1), loop over all shapes
    for (int i = (M-1); i >= 0; --i) {
      
      bool improve = true;
      
      // Possible problem when i=M
      bool test = false;
      if (i < M) test = (shape[i] < shape[i+1]-1);
      
      // Increase i-th shape as long as improvement and not larger than next shape
      while( improve && (i==(M-1) || test) ) {
        
        // Clone vector to avoid changing original one
        IntegerVector new_shape = clone(shape);
        
        // Increase i-th shape by 1
        new_shape[i] = new_shape[i] + 1;
        
        fit = spliceEM_splicefit_raw(pi, theta, new_shape, beta, gamma, 
                                     lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                     trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
        
        new_loglikelihood = as<double>(fit["loglikelihood"]);

        if (new_loglikelihood > loglikelihood + eps) {
          loglikelihood = new_loglikelihood;
          pi = as<double>(fit["pi"]);
          shape = as<IntegerVector>(fit["shape"]);
          theta = as<double>(fit["theta"]);
          beta = as<NumericVector>(fit["beta"]);
          alpha = as<NumericVector>(fit["alpha"]);
          gamma = as<double>(fit["gamma"]);
          
          // number of shape combinations might have changed after EM algorithm if beta_tol > 0
          M = shape.size();
          i = std::min(i, M-1);
          
        } else {
          improve = false;
        }
        
        test = false;
        if (i < M) test = (shape[i] < shape[i+1]-1);
        
      }
      
    }
    
    // Try decreasing the shapes (step 2), loop over all shapes
    for (int i = 0; i<M; i++) {
      
      bool improve = true;
      
      // Possible problem when i=0 or i=M
      bool test = false;
      if (i < M) test = (shape[i] > shape[i-1]+1);
      
      bool test2 = false;
      if (i < M) test2 = (shape[i] > 1);
      
      // Increase i-th shape as long as improvement and larger than 1 and not smaller than previous shape
      while( improve && (i==0 || test) && test2 ) {
        
        // Clone vector to avoid changing original one
        IntegerVector new_shape = clone(shape);
        
        // Decrease i-th shape by 1
        new_shape[i] = new_shape[i] - 1;
        
        fit = spliceEM_splicefit_raw(pi, theta, new_shape, beta, gamma, 
                                     lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                     trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
        
        new_loglikelihood = as<double>(fit["loglikelihood"]);

        if (new_loglikelihood > loglikelihood + eps) {
          loglikelihood = new_loglikelihood;
          pi = as<double>(fit["pi"]);
          theta = as<double>(fit["theta"]);
          shape = as<IntegerVector>(fit["shape"]);
          beta = as<NumericVector>(fit["beta"]);
          alpha = as<NumericVector>(fit["alpha"]);
          gamma = as<double>(fit["gamma"]);
          
          // number of shape combinations might have changed after EM algorithm if beta_tol > 0
          M = shape.size();
          i = std::min(i, M-1);
          
        } else {
          improve = false;
        }
        
        test = false;
        if (i < M) test = (shape[i] > shape[i-1]+1);
        
        test2 = false;
        if (i < M) test2 = (shape[i] > 1);
        
      }
    }
    
    after_loglikelihood = loglikelihood;
    ++iteration;
    
  }  
  
  // Compute ICs
  int df = 2*alpha.size()-1+1+1+1;
  double aic = -2 * loglikelihood + 2.0 * df;
  int n = lower1.size() + lower2.size() + lower3.size() + lower4.size() + lower5.size();
  double bic = -2 * loglikelihood + std::log(n*1.0) * df;
  
  
  // Make output list
  List output = Rcpp::List::create(Rcpp::Named("pi") = pi,
                                   Rcpp::Named("alpha") = alpha,
                                   Rcpp::Named("beta") = beta,
                                   Rcpp::Named("shape") = shape,
                                   Rcpp::Named("theta") = theta,
                                   Rcpp::Named("gamma") = gamma,
                                   Rcpp::Named("loglikelihood") = loglikelihood,
                                   Rcpp::Named("AIC") = aic,
                                   Rcpp::Named("BIC") = bic);
  
  return output;
  
}




// Shape adjustments (export to R)
// [[Rcpp::export]]
List spliceEM_shape_adj_Rexport(const double pi, const double theta, const IntegerVector shape, const NumericVector beta, const double gamma, 
                                const NumericVector lower1, const NumericVector lower2, const NumericVector lower3, const NumericVector lower4, const NumericVector lower5, 
                                const NumericVector upper3, const NumericVector upper4, const NumericVector upper5, 
                                const double trunclower, const double tsplice, const double truncupper,
                                const double eps, const double beta_tol, const double maxiter) {
  
  // Call shape adjustment function
  return spliceEM_shape_adj(pi, theta, shape, beta, gamma, 
                            lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                            trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
}


// Shape reduction (4.3 in Verbelen et al. 2015)
// [[Rcpp::export]]
List spliceEM_shape_red(const double pi_in, const double theta_in, const IntegerVector shape_in, const NumericVector beta_in, const double gamma_in, 
                        const NumericVector lower1, const NumericVector lower2, const NumericVector lower3, const NumericVector lower4, const NumericVector lower5, 
                        const NumericVector upper3, const NumericVector upper4, const NumericVector upper5, 
                        const double trunclower, const double tsplice, const double truncupper, const String criterium, bool improve,
                        const double eps, const double beta_tol, const double maxiter, const bool adj) {
  
  double pi = pi_in;
  double theta = theta_in;
  // Clone vectors to avoid changing it in R (and C++)
  IntegerVector shape = clone(shape_in);
  NumericVector beta = clone(beta_in);
  double gamma = gamma_in;
  
  int n = lower1.size() + lower2.size() + lower3.size() + lower4.size() + lower5.size();
  
  List fit;
  
  // Fit model with initial M
  if (adj) {
    // Adjust shapes
    fit = spliceEM_shape_adj(pi, theta, shape, beta, gamma, 
                             lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                             trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
    
  } else {
    // Do not adjust shapes => apply EM directly
    fit = spliceEM_splicefit_raw(pi, theta, shape, beta, gamma, 
                                 lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                 trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
    
  }

  double loglikelihood = as<double>(fit["loglikelihood"]);
  double IC = as<double>(fit[criterium]);
  pi = as<double>(fit["pi"]);
  theta = as<double>(fit["theta"]);
  shape = as<IntegerVector>(fit["shape"]);
  beta = as<NumericVector>(fit["beta"]);
  NumericVector alpha = as<NumericVector>(fit["alpha"]);
  gamma = as<double>(fit["gamma"]);

  double new_IC;
  
  // Try to improve IC by reducing M
  while( improve && shape.size()>1 ) {

    // Remove shape with smallest beta
    IntegerVector new_shape = shape[beta != min(beta)];
    NumericVector new_beta = beta[beta != min(beta)];
    // Ensure that new beta's sum to 1
    new_beta = new_beta / sum(new_beta);
    
    if (adj) {
      //  Adjust shapes
      fit = spliceEM_shape_adj(pi, theta, new_shape, new_beta, gamma, 
                               lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                               trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
      
    } else {
      // Do not adjust shapes => apply EM directly
      fit = spliceEM_splicefit_raw(pi, theta, new_shape, new_beta, gamma, 
                                   lower1, lower2, lower3, lower4, lower5, upper3, upper4, upper5,
                                   trunclower, tsplice, truncupper, eps, beta_tol, maxiter);
      
    }
    new_IC = as<double>(fit[criterium]);
    
    if (new_IC < IC) {
      IC = new_IC;
      loglikelihood = as<double>(fit["loglikelihood"]);
      pi = as<double>(fit["pi"]);
      theta = as<double>(fit["theta"]);
      shape = as<IntegerVector>(fit["shape"]);
      beta = as<NumericVector>(fit["beta"]);
      alpha = as<NumericVector>(fit["alpha"]);
      gamma = as<double>(fit["gamma"]);
      
    } else {
      improve = false;
    }
  }
  
  // Compute ICs
  int M = alpha.size();
  int df = 2*M-1+1+1+1;
  double aic = -2 * loglikelihood + 2.0 * df;
  double bic = -2 * loglikelihood + std::log(n*1.0) * df;

  
  // Make output list
  List output = Rcpp::List::create(Rcpp::Named("M") = M,
                                   Rcpp::Named("pi") = pi,
                                   Rcpp::Named("alpha") = alpha,
                                   Rcpp::Named("beta") = beta,
                                   Rcpp::Named("shape") = shape,
                                   Rcpp::Named("theta") = theta,
                                   Rcpp::Named("gamma") = gamma,
                                   Rcpp::Named("loglikelihood") = loglikelihood,
                                   Rcpp::Named("AIC") = aic,
                                   Rcpp::Named("BIC") = bic);
  
  return output;
  
}