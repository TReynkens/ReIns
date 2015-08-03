\name{NEWS}
\title{News for package "ReIns"}

\section{Changes in version 0.3-00 (2015-08-03)}{

  \itemize{
    \item Changes in \code{R} code:
      \itemize{
        \item Include GPD: \code{dgpd}, \code{pgpd}, \code{qgpd} and \code{rgpd}.
        \item Include truncated GPD: \code{dtgpd}, \code{ptgpd}, \code{qtgpd} and \code{rtgpd}.
        \item Added functions for splicing of mixed Erlang and EVT distributions.
        \item Added functions for excess-loss of spliced distributions.
        \item Added \code{ES} function to compute expected shortfall.
      }
         
    \item Documentation changes:
      \itemize{
        \item Added documentation for GPD functions.
        \item Added documentation for new splicing functions (including excess-loss).
        \item Added documentation of \code{ES}.
      }
    \item Misscelaneous changes:
      \itemize{
        \item Updated title to have correct upper case format.
        \item Import default packages as required for R 3.3.0.
        \item Add GitHub website as URL for package.
        \item Add bug report page. 
        \item Added Roel Verbelen as author.
      }
  }
}

\section{Changes in version 0.2-01 (2015-06-02)}{
  
  \itemize{
    
    \item Changes in \code{R} code:
       \itemize{
       \item Renamed \code{ExcessMLE}, \code{ProbMLE}, etc. to \code{ExcessGPD}, \code{probGPD}, etc.
      \item Fixed bugs in \code{cEPD}.
      \item Added check in \code{cEPD} that estimates for \eqn{\kappa_1} are larger than \eqn{max(-1,-1/beta.x)} and that \eqn{\gamma} is strictly positive.
      \item Renamed \code{rho1} and \code{beta.x} to \code{rho} and \code{beta} in \code{cEPD}, \code{ciEPD}, \code{cProbEPD}
      and \code{ciProbEPD}.
      \item Removed lines in \code{MeanExcess}.
      \item Added functions \code{cReturn}, \code{ciReturn}, \code{cReturnGH}, \code{cReturnMOM}, \code{cReturnGPD}, \code{cReturnEPD} and \code{ciReturnEPD}.
      \item Added exponential QQplot and log-normal QQplot for right censored data and for interval censored data (\code{cExpQQ}, \code{ciExpQQ}, \code{cLognormalQQ} and \code{ciLognormalQQ}).
      \item Added warning in \code{trciProb} when \eqn{q} is larger than the maximum of \eqn{Z}.
      \item Changed order of argument \code{censored} in all censored probability estimators.
     }
      
      
      
      \item Documentation changes:
             \itemize{
       \item Updated documentation for renaming \code{ExcessMLE}, \code{ProbMLE}, etc. to \code{ExcessGPD}, \code{probGPD}, etc.
         \item Added documentation for \code{cReturn}, \code{ciReturn}, \code{cReturnGH}, \code{cReturnMOM}, \code{cReturnGPD}, \code{cReturnEPD} and \code{ciReturnEPD}.
         \item Added documentation for \code{cExpQQ}, \code{ciExpQQ}, \code{cLognormalQQ} and \code{ciLognormalQQ}.
         \item Updated documentation for changes of argument \code{censored} in all censored probability estimators.
         \item Updated documentation of \code{ExpQQ}, \code{ParetoQQ} and \code{LognormalQQ}.
     }
      
  }
}

\section{Changes in version 0.2-00 (2015-05-20)}{

  \itemize{
  
      \item Changes in \code{R} code:
            \itemize{
            \item Multiple values for \code{rho} can be provided in \code{EPD}.
            \item Fixed bug in \code{pEPD}: it computed the survival function instead of the CDF. 
            \item Fixed bug in \code{ProbEPD} and \code{ReturnEPD}: division and multiplication by the survival function were 
                  switched.
            \item \code{KaplanMeier} can now handle vector input for \code{x} and calls \code{survfit} of \pkg{survival} for 
                  the actual computations if \code{length(X)}>=1000.
            \item The function \code{KaplanMeier.k} is now removed.
            \item New function \code{tr.prob.pos} to estimate exceedance probabilities for truncated data.
            \item Changed \code{gamma} to \code{gamma1} for estimators for censored data.
            \item Added Turnbull estimator for CDF of interval censored data (\code{Turnbull}).
            \item New function \code{tr.cHill} for (right) censored truncated data. Corresponding exceedance probabilities 
                  can be computed using \code{tr.cProb}.
            \item New functions \code{ciHill} for interval censored data and \code{tr.ciHill} for interval censored truncated 
                  data.  Corresponding exceedance probabilities can be computed using using \code{ciProb} and 
                  \code{tr.ciProb}, respectively.
            \item Replaced \code{delta} by \code{kappa} in \code{pEPD} and \code{dEPD} to avoid confusion with indicator for 
                  non-censoring. 
            \item New function \code{ciEPD} to compute the EPD estimator adapted for interval censored data. Corresponding 
                  small exceedance probabilities are computed using \code{ciProbEPD}.
            \item Removed \code{rhoEst} from Namespace.
            \item Added \code{na.rm=TRUE} to \code{min} and \code{max} in \code{checkInput}.
            \item Replaced censored Hill by ordinary Hill in \code{cgenQQ}.
            \item Fixed bug in \code{cGPDmle}: \code{gamma1} was called in output from \code{GPDmle} instead of \code{gamma}.
            \item Fixed bug in \code{ProbGH} and \code{QuantGH}: \code{gamma1} was called in output from \code{Hill} instead 
                  of \code{gamma}. 
            \item Fixed bug in \code{cProbGH} and \code{cProb}: \eqn{(km-1)} should have been \eqn{(1-km)}.
            \item Fixed bug in \code{cProbMLE} and \code{cQuantMLE}: \code{a} is equal to \code{sigma1[K]/pk}, not \code{sigma
                  1[K]^2/pk}.
            \item Added \code{Prob}, \code{Quant}, \code{Quant.2oQV} and \code{Return} as different names for functions
                  \code{Weissman.p}, \code{Weissman.q}, \code{Weissman.q.2oQV} and \code{Weissman.r}.
            \item New functions to compute premiums of excess-loss insurances: \code{Excess}, \code{ExcessEPD} and \code{ExcessMLE}.
          }
          
          
      \item Changes in datasets:
            \itemize{
            \item Fire insurance dataset from the Netherlands: \code{NLfire}.
            }
    
    
     \item Documentation changes:
       \itemize{
       \item Included documentation for \code{Weissman.r} in the documentation for \code{Weissman.p}.
       \item Added documentation for the Pareto and truncated Pareto distribution.
       \item Added documentation for \code{Scale} and \code{Scale.2o}.
       \item Added documentation for \code{EPD}, \code{EPDfit}, \code{ProbEPD} and \code{ReturnEPD}.
       \item Added documentation for \code{KaplanMeier} and \code{Turnbull}.
       \item Added documentation for \code{ProbGH}, \code{ProbMLE}, \code{ProbMOM}, \code{QuantGH}, \code{QuantMLE} and     
              \code{QuantMOM}.
       \item Added documentation for \code{ParetoQQ_der} and \code{LognormalQQ_der}.
       \item Updated documentation for \code{ExpQQ} and \code{MeanExcess}.
       \item Added documentation for \code{ScaleReg}, \code{ProbReg} and \code{QuantReg}.
       \item Added documentation for \code{GPDresiduals}.
       \item Added documentation for \code{cEPD}, \code{ciEPD}, \code{cProbEPD} and
             \code{ciProbEPD}.
        \item Added documentation for \code{cHill}, \code{ciHill}, \code{cProb}, \code{cQuant} and
             \code{ciProb}.
        \item Added documentation for \code{cgenHill}, \code{cParetoQQ} and  \code{cgenQQ}.
        \item Added documentation for \code{cMoment} and\code{cGPDmle}.
        \item Added documentation for \code{cProbGH}, \code{cProbMOM}, \code{cProbMLE}, \code{cQuantGH},\code{cQuantMOM} and 
              \code{cQuantMLE}.
        \item Added documentation for \code{trHill}, \code{trDT}, \code{trEndpoint}, \code{trProb} and \code{trQuant}.
        \item Updated documentation of \code{Weissman.p}, \code{Weissman.q}, \code{Weissman.q.2oQV} and \code{Weissman.r}
        \item Added documentation for \code{trcHill}, \code{trciHill}, \code{trcProb} and \code{trciProb}.
        \item Added documentation for \code{Excess}, \code{ExcessEPD} and \code{ExcessMLE}.
       }
  }
}

\section{Changes in version 0.1-01 (2015-04-30)}{

  \itemize{
  
      \item Changes in \code{R} code:
            \itemize{
            \item New function \code{EPD} which computes the EPD estimator of Beirlant et al. (2009).
            The corresponding estimator of the return levels is implemented in \code{returnEPD}.
            \item New function \code{rhoEst} which computes the estimator for \eqn{\rho} of Fraga Alves et al. (2003).
            \item New functions \code{ScaleReg}, \code{ScaleReg_surv} and \code{ScaleReg_survi} to compute scale
                  estimates in a regression context with \eqn{\gamma} constant. The functions \code{ProbReg} and \code{QuantReg
                } compute the corresponding estimates of the small exceedance probabilities and large quantiles, respectively.
            \item The function \code{checkInput} now also (optionally) checks the input arguments \code{scale} and \code{DT}.
            \item New functions \code{checkProb} and \code{checkCensored} to check if an argument is a valid probability or a                    valid logical vector of the needed length, respectively.
            \item New functions \code{cProbGH}, \code{cQuantGH}, \code{cProbMOM}, \code{cQuantMOM}, \code{cProbMLE} and                             \code{cQuantMLE} to compute small exceedance probabilities and large quantiles using the censored
                      generalised Hill, moment and POT (MLE) estimators, respectively.
            \item Added input checking (using \code{checkProb}) for the input argument \code{p} in all functions that estimate 
                  large quantiles.     
            \item All estimators for small exceedance probabilities now return \code{NA} when an estimated probability is not a                   valid probability (smaller than 0 or larger than 1 ).
            \item New function \code{cPOT} to compute censored POT estimates.
            \item New function \code{KaplanMeier} to compute the Kaplan-Meier estimator for the CDF.
            \item Fixed bugs with Kaplan-Meier estimator in \code{cParetoQQ}, \code{cgenQQ}, \code{cProbGH}, \code{cQuantGH}, 
                  \code{cProbMOM}, \code{cQuantMOM}, \code{cProbMLE} and \code{cQuantMLE}.
            \item New functions \code{cEPD} and \code{cProbEPD} to compute the EPD estimator and corresponding small exceedance                   probabilities for censored data, (Bardoutsos et al., 2015).
            }
  
      
 
      \item Changes in testing:
         \itemize{
            \item Included test for \code{p.hat.fun}.
            }

  
  }
  
}


