#ifndef SpliceEM_Mstep_H
#define SpliceEM_Mstep_H


#include "SpliceEM_aux.h"

// Compute theta in M-step of ME
double spliceEM_theta(const double ltheta, const double E1_ME, const double E3_ME, const double E5_ME,
                      const double n1_h, const NumericVector beta, const NumericVector shape, const double trunclower, const double tsplice);

// Compute gamma in M-step of Pareto
double spliceEM_Mstep_Pareto(const double lgamma, const double H, const double tsplice, const double truncupper);

#endif