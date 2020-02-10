\name{NEWS}
\title{News for package "ReIns"}
\encoding{UTF-8}


\section{Changes in version 1.0.9 (2020-02-10)}{

   \subsection{Changes in \code{R} code:}{
     \itemize{
      \item Fix bug in \code{cgenHill}.
      \item Fix bug in \code{cQuantGH}.
    }
   }
}


\section{Changes in version 1.0.8 (2018-09-18)}{

   \subsection{Changes in \code{C++} code:}{
     \itemize{
      \item Load \code{cfloat} to define \code{DBL_MAX} when strict R headers are set (flag \code{-DSTRICT_R_HEADERS} in MAKEVARS) as this will be done by default in \pkg{Rcpp} in September 2019. See \url{https://github.com/RcppCore/Rcpp/issues/898} for more details. This fixes \url{https://github.com/TReynkens/ReIns/issues/2}.
    }
   }
   
   \subsection{Miscellaneous changes:}{
     \itemize{
      \item Update URL for SOA data.
    }
  }      
}

\section{Changes in version 1.0.7 (2018-04-24)}{


  \subsection{Changes in \code{R} code:}{
    \itemize{
     \item Fix bug with random samples of Mixed Erlang when only one shape is present.
     \item Ensure that \code{rSplice} also works for the ME-GPD model.
    }
  }
  
   \subsection{Miscellaneous changes:}{
     \itemize{
      \item Add ORCID identifiers.
    }
 }      

}

\section{Changes in version 1.0.6 (2017-12-11)}{

  \subsection{Changes in documentation:}{
    \itemize{
      \item Update Reynkens et al. (2017) reference to final version.
    }
  }      
   \subsection{Miscellaneous changes:}{
     \itemize{
      \item Correct link to Wiley webpage for Albrecher et al. (2017). 
    }
 }      

}

\section{Changes in version 1.0.5 (2017-09-07)}{

  \subsection{Changes in \code{C++} code:}{
    \itemize{
      \item Change calls to \code{C++} functions as required by \pkg{Rcpp} version 0.12.12.
    }
  } 

  \subsection{Changes in documentation:}{
    \itemize{
      \item Update Reynkens et al. (2017) reference.
    }
  }      
        
  \subsection{Miscellaneous changes:}{
     \itemize{
      \item Add link to Unil webpage for Albrecher et al. (2017) in DESCRIPTION.
      \item Enable byte-compilation of package code.
      \item Set minimal required \pkg{Rcpp} version to 0.12.12.
    }
 }
}


\section{Changes in version 1.0.4 (2017-06-10)}{


  \subsection{Changes in \code{R} code:}{
    \itemize{
       \item Remove unnecessary capitals in plot labels.
       \item Make sure that both \code{R} and \code{L} have length 1 in \code{ExcessEPD}, \code{ExcessGPD}
    and \code{ExcessPareto}.
    }
  }
  
  \subsection{Changes in documentation:}{
    \itemize{
      \item Include default value for \code{L} in documentation of \code{ExcessSplice}.
      \item Update Reynkens et al. (2017) reference.
      \item Update Beirlant et al. (2017) reference.
    }
  }      
  
  \subsection{Miscellaneous changes:}{
    \itemize{
      \item Add reference to website of Albrecher et al. (2017) in DESCRIPTION.
    }
  }
}


\section{Changes in version 1.0.3 (2017-04-18)}{

 \itemize{
    \item Add estimators from Beirlant et al. (2017) that can handle truncated distributions with \eqn{\gamma>-0.5}.
    \item \code{C++}: use output from \code{package_native_routine_registration_skeleton()} (\pkg{tools}) to 
                      register symbols and disable symbol search.
 }
}

\section{Changes in version 1.0.2 (2017-02-15)}{

 \itemize{
    \item \code{C++}: register native symbols and disable symbol search (as requested by NOTE in R CMD check).
 }
}

\section{Changes in version 1.0.1 (2016-12-26)}{

 \itemize{
    \item Fix overloaded use of logarithm with integers in \code{C++} to avoid errors on Solaris (as requested by CRAN).
    \item Change reference to Albrecher et al. (2017) as requested by CRAN.
    \item Change vignette index entry (shown on CRAN).
 }
}

\section{Changes in version 1.0.0 (2016-12-23)}{

 \itemize{
    \item Release on CRAN. The old version of \code{NEWS.Rd} can be found on GitHub: \url{https://github.com/TReynkens/ReIns/blob/master/inst/NEWS.Rd}.
 }
}


