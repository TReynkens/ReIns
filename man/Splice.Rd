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
Output list from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitciPareto}} or \code{\link{SpliceFitGPD}}.
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
See Section 3.3 of Albrecher et al. (2016) for details.
}
\value{
Vector of densities \eqn{f(x)} for \code{dSplice}, vector of probabilities \eqn{F(x)=P(X\le x)} for \code{pSplice},  vector of quantiles \eqn{Q(p)} for \code{qSplice} and a random sample of length \eqn{n} for \code{rSplice}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang CDF, PDF and quantiles.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{VaR}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitciPareto}}, \code{\link{SpliceFitGPD}},
\code{\link{SpliceECDF}}, \code{\link{SpliceLL}}, \code{\link{SplicePP}}
}

\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)



x <- seq(0,20,0.01)

# Plot of splicing CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of splicing PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")



p <- seq(0,1,0.01)

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

