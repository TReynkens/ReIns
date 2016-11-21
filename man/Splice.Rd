\name{Splice}
\alias{dSplice}
\alias{pSplice}
\alias{qSplice}
\alias{rSplice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spliced distribution
}
\description{
Density, distribution function, quantile function and random generation for the fitted spliced distribution.
}
\usage{
dSplice(x, splicefit, log = FALSE)

pSplice(x, splicefit, lower.tail = TRUE, log.p = FALSE)

qSplice(p, splicefit, lower.tail = TRUE, log.p = FALSE)

rSplice(n, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to evaluate the CDF or PDF in.
}
  \item{p}{
Vector of probabilities.
}
  \item{n}{
Number of observations.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
 \item{log}{
  Logical indicating if the densities are given as \eqn{\log(f)}, default is \code{FALSE}.
 }
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as \eqn{\log(p)}, default is \code{FALSE}.
 }
}
\details{ 
See Reynkens et al. (2016) and Section 4.3 in Albrecher et al. (2017) for details.
}
\value{
\code{dSplice} gives the density function evaluated in \eqn{x}, \code{pSplice} the CDF evaluated in \eqn{x} and \code{qSplice} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \eqn{x} or \eqn{p}.

\code{rSplice} returns a random sample of length \eqn{n}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2016). "Modeling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang PDF, CDF and quantiles.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{VaR}}, \code{\link{SpliceFit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceFitGPD}},
\code{\link{SpliceECDF}}, \code{\link{SpliceLL}}, \code{\link{SplicePP}}
}

\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)



x <- seq(0, 20, 0.01)

# Plot of spliced CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of spliced PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")



p <- seq(0, 1, 0.01)

# Plot of splicing quantiles
plot(p, qSplice(p, splicefit), type="l", xlab="p", ylab="Q(p)")

# Plot of VaR
plot(p, VaR(p, splicefit), type="l", xlab="p", ylab=bquote(VaR[p]))



# Random sample from spliced distribution
x <- rSplice(1000, splicefit)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

