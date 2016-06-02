\name{GPD}
\alias{dgpd}
\alias{pgpd}
\alias{qgpd}
\alias{rgpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The generalised Pareto distribution
}
\description{
Density, distribution function, quantile function and random generation for the Generalised Pareto Distribution (GPD).
}
\usage{
dgpd(x, gamma, mu = 0, sigma, log = FALSE)
pgpd(x, gamma, mu = 0, sigma, lower.tail = TRUE, log.p = FALSE)
qgpd(p, gamma, mu = 0, sigma, lower.tail = TRUE, log.p = FALSE)
rgpd(n, gamma, mu = 0, sigma)
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
 \item{log}{
  Logical indicating if the densities are given as log(f), default is \code{FALSE}.
 }
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as \eqn{log(p)}, default is \code{FALSE}.
 }
}
\details{
The Cumulative Distribution Function (CDF) of the GPD for \eqn{\gamma != 0} is equal to
\eqn{F(x) = 1-(1+\gamma (x-\mu)/\sigma)^{-1/\gamma}} for all \eqn{x \ge \mu} and \eqn{F(x)=0} otherwise.
When \eqn{\gamma=0}, the CDF is given by 
\eqn{F(x) = 1-\exp((x-\mu)/\sigma)} for all \eqn{x \ge \mu} and \eqn{F(x)=0} otherwise.
}
\value{
\code{dgpd} gives the density function evaluated in \eqn{x}, \code{pgpd} the CDF evaluated in \code{x} and \code{qgpd} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rgpd} returns a random sample of length \eqn{n}.
}

\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ptgpd]{tGPD}}, \code{\link[=ppareto]{Pareto}}, \code{\link[=pepd]{EPD}}, \code{\link[stats]{Distributions}}
}


\examples{
# Plot of the PDF
x = seq(0,10,0.01)
plot(x, dgpd(x, gamma=1/2, sigma=5), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(0,10,0.01)
plot(x, pgpd(x, gamma=1/2, sigma=5), xlab="x", ylab="CDF", type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}