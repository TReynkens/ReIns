\name{NEWS}
\title{News for package "ReIns"}
\encoding{UTF-8}

\section{Changes in version 1.0.5 (2017-yy-xx)}{

 \itemize{
    \item Add link to Unil webpage for Albrecher et al. (2017) in DESCRIPTION.
 }
}


\section{Changes in version 1.0.4 (2017-06-10)}{

 \itemize{
    \item Remove unnecessary capitals in plot labels.
    \item Update Reynkens et al. (2017) reference.
    \item Update Beirlant et al. (2017) reference.
    \item Include default value for \code{L} in documentation of \code{ExcessSplice}.
    \item Make sure that both \code{R} and \code{L} have length 1 in \code{ExcessEPD}, \code{ExcessGPD}
    and \code{ExcessPareto}.
    \item Add reference to website of Albrecher et al. (2017) in DESCRIPTION.
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


