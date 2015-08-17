\name{SpliceCDF}
\alias{SpliceCDF}
\alias{SplicePDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CDF and PDF of splicing fit
}
\description{
Compute the Cumulative Distribution Function (CDF) and Probability Density Function (PDF) of the fitted splicing distribution.
}
\usage{
SpliceCDF(x, splicefit, lower.tail = TRUE, log.p = FALSE)

SplicePDF(x, splicefit, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to evaluate the CDF or PDF in.
}
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}} or \code{\link{SpliceFitGPD}}.
}
 \item{log}{
  Logical indicating if the densities are given as log(f), default is \code{FALSE}.
 }
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as log(p), default is \code{FALSE}.
 }
}
\details{ 
See Section 3.3 of Albrecher et al. (2015) for details.
}
\value{
Vector of probabilities \eqn{F(x)=P(X\le x)} for \code{SpliceCDF} and vector of densities \eqn{f(x)} for \code{SplicePDF}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang CDF and PDF.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceQuant}}, \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}},
\code{\link{SpliceECDF}}, \code{\link{SpliceLL}}, \code{\link{SplicePP}}
}

\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.6)

x <- seq(0,20,0.01)
# Plot of splicing CDF
plot(x, SpliceCDF(x, splicefit), type="l", xlab="x", ylab="F(x)")
# Plot of splicing PDF
plot(x, SplicePDF(x, splicefit), type="l", xlab="x", ylab="f(x)")

# Fitted survival function and empirical survival function 
SpliceECDF(x, X, splicefit)

# Log-log plot with empirical survival function and fitted survival function
SpliceLL(x, X, splicefit)

# PP-plot of empirical survival function and fitted survival function
SplicePP(x, X, splicefit)

# PP-plot of empirical survival function and 
# fitted survival function with log-scales
SplicePP(x, X, splicefit,log=TRUE)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

