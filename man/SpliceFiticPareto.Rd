\name{SpliceFiticPareto}
\alias{SpliceFiticPareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing of mixed Erlang and Pareto for interval censored data
}

\description{
Fit spliced distribution of a mixed Erlang distribution and a Pareto distribution adapted for interval censoring and truncation.
}
\usage{
SpliceFiticPareto(L, U, censored, tsplice, M = 3, s = 1:10, trunclower = 0, 
                  truncupper = Inf, ncores = NULL, criterium = c("BIC", "AIC"), 
                  reduceM = TRUE, eps = 10^(-3), beta_tol = 10^(-5), maxiter = Inf, 
                  cpp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
Vector of length \eqn{n} with the lower boundaries of the intervals for interval censored data or the observed data for right censored data.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundaries of the intervals.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{tsplice}{
The splicing point \eqn{t}.
}

  \item{M}{
Initial number of Erlang mixtures, default is 3. This number can change when determining
an optimal mixed Erlang fit using an information criterion.
}
  \item{s}{
Vector of spread factors for the EM algorithm, default is \code{1:10}. We loop over these factors when determining
an optimal mixed Erlang fit using an information criterion, see Verbelen et al. (2016).
}
  \item{trunclower}{
Lower truncation point. Default is 0.
}
  \item{truncupper}{
Upper truncation point. Default is \code{Inf} (no upper truncation).
}
  \item{ncores}{
Number of cores to use when determining an optimal mixed Erlang fit using an information criterion.
When \code{NULL} (default), \code{max(nc-1,1)} cores are used where \code{nc} is the number of cores as determined by \code{\link{detectCores}}.
}
\item{criterium}{
Information criterion used to select the number of components of the ME fit and \code{s}. One of \code{"AIC"} and \code{"BIC"} (default).
}
\item{reduceM}{
Logical indicating if M should be reduced based on the information criterion, default is \code{TRUE}.
}
\item{eps}{
Covergence threshold used in the EM algorithm. Default is \code{10^(-3)}.
}
\item{beta_tol}{
Threshold for the mixing weights below which the corresponding shape parameter vector is considered neglectable (ME part). Default is \code{10^(-5)}.
}
 \item{maxiter}{
Maximum number of iterations in a single EM algorithm execution. Default is \code{Inf} meaning no maximum number of iterations.
}
 \item{cpp}{
Use \code{C++} implementation (\code{cpp=TRUE}) or \code{R} implementation (\code{cpp=FALSE}) of the algorithm. Default is \code{FALSE} meaning the plain \code{R} implementation is used.
}
}
\details{
See Reynkens et al. (2017), Section 4.3.2 of Albrecher et al. (2017) and Verbelen et al. (2015) for details. The code follows the notation of the latter. Initial values follow from Verbelen et al. (2016).

Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and left censored data should be entered as \code{L=trunclower} and \code{U=u}. 
}
\value{
A \code{\link{SpliceFit}} object.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758.

Verbelen, R., Antonio, K. and Claeskens, G. (2016). "Multivariate Mixtures
of Erlangs for Density Estimation Under Censoring." \emph{Lifetime Data Analysis}, 22, 429--455.
}
\author{
Tom Reynkens based on \code{R} code from Roel Verbelen for fitting the mixed Erlang distribution (without splicing).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFitPareto}}, \code{\link{SpliceFitGPD}}, \code{\link[=pSplice]{Splice}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X,Y)

# Censoring indicator
censored <- (X>Y)

# Right boundary
U <- Z
U[censored] <- Inf

# Splice ME and Pareto
splicefit <- SpliceFiticPareto(L=Z, U=U, censored=censored, tsplice=quantile(Z,0.9))



x <- seq(0,20,0.1)

# Plot of spliced CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of spliced PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")


# Fitted survival function and Turnbull survival function 
SpliceTB(x, L=Z, U=U, censored=censored, splicefit=splicefit)


# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB(x, L=Z, U=U, censored=censored, splicefit=splicefit)


# PP-plot of Turnbull survival function and fitted survival function
SplicePP_TB(L=Z, U=U, censored=censored, splicefit=splicefit)

# PP-plot of Turnbull survival function and 
# fitted survival function with log-scales
SplicePP_TB(L=Z, U=U, censored=censored, splicefit=splicefit, log=TRUE)

# QQ-plot using Turnbull survival function and fitted survival function
SpliceQQ_TB(L=Z, U=U, censored=censored, splicefit=splicefit)
}
}


