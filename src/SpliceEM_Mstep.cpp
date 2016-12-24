#include "SpliceEM_Mstep.h"

/////////////////////////////////////////
// M-step ME

// Function to compute root of in spliceEM_theta
double f_theta(const double ltheta, const double D, const double br_sum, const NumericVector beta, const IntegerVector shape,
               const double trunclower, const double tsplice) {

  double theta = exp(ltheta);
  
  int M = beta.size();

  double f = theta - D;

  double Ctheta;
  double Ctheta2;


  // Avoid NaN
  if (R_FINITE(tsplice)) {

    for(int j = 0; j < M; j++) {

      // Take logs for numerical stability
      Ctheta = shape[j] * log(trunclower) - trunclower / theta - (shape[j]-1) * log(theta) - R::lgammafn(shape[j]) - log(pGamma(tsplice, shape[j], theta) - pGamma(trunclower, shape[j], theta));

      Ctheta2 = shape[j] * log(tsplice) - tsplice / theta - (shape[j]-1) * log(theta) - R::lgammafn(shape[j]) - log(pGamma(tsplice, shape[j], theta) - pGamma(trunclower, shape[j], theta));

      Ctheta = exp(Ctheta) - exp(Ctheta2);

      // // Without log
      // Ctheta = pow(trunclower, shape[j]) * exp(-trunclower/theta);
      // Ctheta -= pow(tsplice, shape[j]) * exp(-tsplice/theta);
      // Ctheta /= pow(theta, shape[j]-1) * R::gammafn(shape[j]) * (pGamma(tsplice, shape[j], theta) - pGamma(trunclower, shape[j], theta));

      f += beta[j] * Ctheta / br_sum;

    }

  } else {

    for(int j = 0; j < M; j++) {
      // Take logs for numerical stability
      Ctheta = shape[j] * log(trunclower) - trunclower / theta - (shape[j]-1) * log(theta) - R::lgammafn(shape[j]) - log(1 - pGamma(trunclower, shape[j], theta));

      Ctheta = exp(Ctheta);

      f += beta[j] * Ctheta / br_sum;
    }
  }

  return f;
}

// Derivative of f_theta
double f_theta_der(const double ltheta, const double br_sum, const NumericVector beta, const IntegerVector shape,
                   const double trunclower, const double tsplice) {

  double theta = exp(ltheta);
  
  int M = beta.size();

  double fd = 1;

  double dCtheta1;
  double dCtheta2;
  double dCtheta3;
  double dCtheta_tmp1;
  double dCtheta_tmp2;
  double FE_tl;
  double FE_t;

  // Avoid NaN
  if (R_FINITE(tsplice)) {

    for(int j = 0; j < M; j++) {

      FE_tl = pGamma(trunclower, shape[j], theta);
      FE_t = pGamma(tsplice, shape[j], theta);

      // Take log for numerical stability
      dCtheta_tmp1 = (shape[j]+1) * log(theta) +  R::lgammafn(shape[j]) + log(FE_t - FE_tl);

      dCtheta1 = exp((shape[j]+1) * log(trunclower) - trunclower / theta - dCtheta_tmp1);
      dCtheta1 -= exp((shape[j]+1) * log(tsplice) - tsplice / theta - dCtheta_tmp1);

      dCtheta_tmp1 = (shape[j]-1) * log(theta) +  R::lgammafn(shape[j]) + log(FE_t - FE_tl);
      dCtheta2 = shape[j] * log(trunclower) - trunclower / theta - 2 * dCtheta_tmp1;
      dCtheta3 = shape[j] * log(tsplice) - tsplice / theta - 2 * dCtheta_tmp1;

      dCtheta_tmp2 = exp(log(shape[j]-1.0) + R::lgammafn(shape[j]) + (shape[j]-2) * log(theta)  + log(FE_t - FE_tl));
      dCtheta_tmp2 += exp(shape[j] * log(trunclower)  - trunclower / theta - 2 * log(theta)) - exp(shape[j] * log(tsplice) - tsplice / theta - 2 * log(theta));

      dCtheta2 = (exp(dCtheta2) - exp(dCtheta3)) * dCtheta_tmp2;


      // // Without log
      // dCtheta1 = pow(trunclower, shape[j]+1) * exp(-trunclower/theta);
      // dCtheta1 -= pow(tsplice, shape[j]+1) * exp(-tsplice/theta);
      // dCtheta1 /= pow(theta, shape[j]+1) * R::gammafn(shape[j]) * (FE_t - FE_tl);
      // 
      // dCtheta2 = pow(trunclower, shape[j]) * exp(-trunclower/theta);
      // dCtheta2 -= pow(tsplice, shape[j]) * exp(-tsplice/theta);
      // dCtheta2 *= (shape[j]-1) * R::gammafn(shape[j]) * pow(theta, shape[j]-2) * (FE_t - FE_tl) + (pow(trunclower, shape[j])  * exp(-trunclower/theta) - pow(tsplice, shape[j])  * exp(-tsplice/theta) ) / pow(theta, 2.0);
      // dCtheta2 /= pow(pow(theta, shape[j]-1) * R::gammafn(shape[j]) * (FE_t - FE_tl), 2.0);

      fd += beta[j] * (dCtheta1 - dCtheta2) / br_sum;

    }

  } else {

    for(int j = 0; j < M; j++) {

      FE_tl = pGamma(trunclower, shape[j], theta);
      // Take log for numerical stability
      dCtheta_tmp1 = (shape[j]+1) * log(theta) +  R::lgammafn(shape[j]) + log(1 - FE_tl);

      dCtheta1 = exp((shape[j]+1) * log(trunclower) - trunclower / theta - dCtheta_tmp1);

      dCtheta_tmp1 = (shape[j]-1) * log(theta) +  R::lgammafn(shape[j]) + log(1 - FE_tl);
      dCtheta2 = shape[j] * log(trunclower) - trunclower / theta - 2 * dCtheta_tmp1;

      dCtheta_tmp2 = exp(log(shape[j]-1.0) + R::lgammafn(shape[j]) + (shape[j]-2) * log(theta)  + log(1 - FE_tl));
      dCtheta_tmp2 += exp(shape[j] * log(trunclower)  - trunclower / theta - 2 * log(theta));

      dCtheta2 = exp(dCtheta2) * dCtheta_tmp2;

      fd += beta[j] * (dCtheta1 - dCtheta2) / br_sum;

    }

  }

  // Extra part in derivative since log(gamma) is considered
  return fd * theta;
}




// Compute log(theta) in M-step of ME
double spliceEM_theta(const double ltheta, const double E1_ME, const double E3_ME, const double E5_ME,
                      const double n1_h, const NumericVector beta, const IntegerVector shape, const double trunclower, const double tsplice) {
  
  int M = beta.size();
  
  // Sum of beta*shape
  double br_sum = 0;
  
  for(int j = 0; j < M; j++) {
    
    br_sum += beta[j]*shape[j];
  }

  double D = (E1_ME + E3_ME + E5_ME) / (n1_h * br_sum);

  // Newton-Raphson approach
  
  // Number of iterations-1
  int iter = 0;
  // Maximal number of iterations
  int max_iter = 100;
  // Numerical tolerance
  double tol = pow(10, -6.0);

  // Initial guess
  double m_old = ltheta;
  double m = ltheta;

  double f;
  double f_der;
  
  // No optimisation needed when no truncation
  if (trunclower==0 && !R_FINITE(tsplice)) {
    m = log(D);

  } else {

    // Optimisation loop
    while ( (std::abs((m-m_old)/m_old) > tol || iter==0) && iter < max_iter) {

      m_old = m;


      f = f_theta(m_old, D, br_sum, beta, shape, trunclower, tsplice);
      f_der = f_theta_der(m_old, br_sum, beta, shape, trunclower, tsplice);

      // Problems when f is infinite or NaN
      if (!R_FINITE(f) || R_IsNaN(f)) {
        f = DBL_MAX;
      }

      // Problems when derivative is zero
      if (std::abs(f_der) < pow(10, -14.0)) break;

      // Problems when derivative is infinite or NaN
      if (!R_FINITE(f_der) || R_IsNaN(f_der)) {
        f_der = DBL_MAX;
      }

      // New estimate
      m = m_old - f / f_der;

      // Increase iteration
      ++iter;
    }

  }

  return m;
}


/////////////////////////////////////////
// M-step Pareto



// Function to compute root of in spliceEM_Mstep_Pareto
double f_Mgamma (const double lgamma, const double H, const double beta_trunc) {
  double gamma = exp(lgamma);
  return gamma - H - log(beta_trunc) / (pow(beta_trunc, 1/gamma) - 1);
}


// Derivative of f_trHill
double f_Mgamma_der (const double lgamma, const double H, const double beta_trunc) {
  double gamma = exp(lgamma);
  return 1 - 1/pow(gamma, 2.0) * pow(log(beta_trunc), 2.0) * pow(beta_trunc, 1/gamma) / pow(pow(beta_trunc, 1/gamma)-1, 2.0);
  // Take logs for numerical stability
  //return 1 - exp(-2*log(gamma) + 2*log(log(beta_trunc)) + 1/gamma*log(beta_trunc) - 2*log(pow(beta_trunc, 1/gamma)-1));
}



// Compute log(gamma) in M-step of Pareto
// gamma is the starting value, H the estimate for gamma without truncation
double spliceEM_Mstep_Pareto(const double lgamma, const double H, const double tsplice, const double truncupper) {
  
  const double beta_trunc = truncupper / tsplice;
  // Newton-Raphson approach
  
  // Number of iterations-1
  int iter = 0;
  // Maximal number of iterations
  int max_iter = 100;
  // Numerical tolerance
  double tol = pow(10, -6.0);

  // Initial guess
  double m_old = lgamma;
  double m = lgamma;
  
  double f;
  double f_der;
  
  // Optimisation loop
  while ( (std::abs((m-m_old)/m_old) > tol || iter==0) && iter < max_iter) {
  
    m_old = m;
    
    f = f_Mgamma(m_old, H, beta_trunc);
    // Extra part in derivative since log(gamma) is considered
    f_der = f_Mgamma_der(m_old, H, beta_trunc) * exp(m_old);
    
    // Problems when f is infinite or NaN
    if (!R_FINITE(f) || R_IsNaN(f)) {
      f = DBL_MAX;
    }
    
    // Problems when derivative is zero
    if (std::abs(f_der) < pow(10, -14.0)) break;
    
    // Problems when derivative is infinite or NaN
    if (!R_FINITE(f_der) || R_IsNaN(f_der)) {
      f_der = DBL_MAX;
    }
    
    // New estimate
    m = m_old -  f / f_der;
    
    // Increase iteration
    ++iter;
  } 

  return m;
  
}


