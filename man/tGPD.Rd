\name{tGPD}
\alias{dtgpd}
\alias{ptgpd}
\alias{qtgpd}
\alias{rtgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The truncated generalised Pareto distribution
}
\description{
Density, distribution function, quantile function and random generation for the truncated Generalised Pareto Distribution (GPD).
}
\usage{
dtgpd(x, gamma, mu = 0, sigma, endpoint = Inf, log = FALSE)
ptgpd(x, gamma, mu = 0, sigma, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
qtgpd(p, gamma, mu = 0, sigma, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
rtgpd(n, gamma, mu = 0, sigma, endpoint = Inf)
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
The \eqn{\gamma} parameter of the GPD, a real number.
}
  \item{mu}{
The \eqn{\mu} parameter of the GPD, a strictly positive number. Default is 0.
}
  \item{sigma}{
The \eqn{\sigma} parameter of the GPD, a strictly positive number.
}
  \item{endpoint}{
Endpoint of the truncated GPD. The default value is \code{Inf} for which the truncated GPD corresponds to the ordinary GPD.
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
The Cumulative Distribution Function (CDF) of the truncated GPD is equal to
\eqn{F_T(x) = F(x) / F(T)} where \eqn{F} is the CDF of the ordinary GPD and \eqn{T} is the endpoint (truncation point) of the truncated GPD.
}
\value{
\code{dtgpd} gives the density function evaluated in \eqn{x}, \code{ptgpd} the CDF evaluated in \eqn{x} and \code{qtgpd} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \eqn{x} or \eqn{p}.

\code{rtgpd} returns a random sample of length \eqn{n}.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link[=ptgpd]{tGPD}}, \code{\link[=ppareto]{Pareto}}, \code{\link[stats]{Distributions}}
}


\examples{
# Plot of the PDF
x <- seq(0, 10, 0.01)
plot(x, dtgpd(x, gamma=1/2, sigma=5, endpoint=8), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x <- seq(0, 10, 0.01)
plot(x, ptgpd(x, gamma=1/2, sigma=5, endpoint=8), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}