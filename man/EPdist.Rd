\name{Extended Pareto}
\alias{depd}
\alias{pepd}
\alias{qepd}
\alias{repd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Extended Pareto Distribution
}
\description{
Density, distribution function, quantile function and random generation for the Extended Pareto Distribution (EPD).
}
\usage{
depd(x, gamma, kappa, tau = -1, log = FALSE)
pepd(x, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE)
qepd(p, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE)
repd(n, gamma, kappa, tau = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
Vector of quantiles.
}
  \item{p}{
Vector of probabilities.
}
  \item{n}{
Number of observations.
}

  \item{gamma}{
The \eqn{\gamma} parameter of the EPD, a strictly positive number.
}
  \item{kappa}{
The \eqn{\kappa} parameter of the EPD. It should be larger than \eqn{\max{-1,1/\tau}}.
}
  \item{tau}{
The \eqn{\tau} parameter of the EPD, a strictly negative number. Default is -1.
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
The Cumulative Distribution Function (CDF) of the EPD for is equal to
\eqn{F(x) = 1-(x(1+\kappa-\kappa x^{\tau})^{-1/\gamma}} for all \eqn{x > 1} and \eqn{F(x)=0} otherwise.

Note that an EPD random variable with \eqn{\tau=-1} and \eqn{\kappa=\gamma/\sigma-1} is GPD distributed with \eqn{\mu=1}, \eqn{\gamma} and \eqn{\sigma}.
}
\value{
\code{depd} gives the density function evaluated in \eqn{x}, \code{pepd} the CDF evaluated in \code{x} and \code{qepd} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{repd} returns a random sample of length \eqn{n}.
}

\references{
Beirlant, J., Joossens, E. and Segers, J. (2009). "Second-Order Refined Peaks-Over-Threshold
Modelling for Heavy-Tailed Distributions." \emph{Journal of Statistical Planning and Inference}, 139, 2800--2815.

}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ppareto]{Pareto}}, \code{\link[=pgpd]{GPD}}, \code{\link[stats]{Distributions}}
}


\examples{
# Plot of the PDF
x <- seq(0, 10, 0.01)
plot(x, depd(x, gamma=1/2, kappa=1, tau=-1), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x <- seq(0, 10, 0.01)
plot(x, pepd(x, gamma=1/2, kappa=1, tau=-1), xlab="x", ylab="CDF", type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
