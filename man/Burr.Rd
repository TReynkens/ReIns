\name{Burr}
\alias{dburr}
\alias{pburr}
\alias{qburr}
\alias{rburr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Burr distribution
}
\description{
Density, distribution function, quantile function and random generation for the Burr distribution (type XII).
}
\usage{
dburr(x, alpha, rho, eta = 1, log = FALSE)
pburr(x, alpha, rho, eta = 1, lower.tail = TRUE, log.p = FALSE)
qburr(p, alpha, rho, eta = 1, lower.tail = TRUE, log.p = FALSE)
rburr(n, alpha, rho, eta = 1)
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
  \item{alpha}{
The \eqn{\alpha} parameter of the Burr distribution, a strictly positive number.
}
  \item{rho}{
The \eqn{\rho} parameter of the Burr distribution, a strictly negative number.
}
  \item{eta}{
The \eqn{\eta} parameter of the Burr distribution, a strictly positive number.
The default value is \code{1}.
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
The Cumulative Distribution Function (CDF) of the Burr distribution is equal to
\eqn{F(x) = 1-((\eta+x^{-\rho*\alpha})/\eta)^{1/\rho}} for all \eqn{x \ge 0} and \eqn{F(x)=0} otherwise. We need that \eqn{\alpha>0}, \eqn{\rho<0} and \eqn{\eta>0}.

Beirlant et al. (2004) uses parameters \eqn{\eta, \tau, \lambda} which correspond to \eqn{\eta},  \eqn{\tau=-\rho*\alpha} and \eqn{\lambda=-1/\rho}.
}
\value{
\code{dburr} gives the density function evaluated in \eqn{x}, \code{pburr} the CDF evaluated in \code{x} and \code{qburr} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rburr} returns a random sample of length \eqn{n}.
}

\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ptburr]{tBurr}},  \code{\link[stats]{Distributions}}
}
\examples{
# Plot of the PDF
x <- seq(0, 10, 0.01)
plot(x, dburr(x, alpha=2, rho=-1), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x <- seq(0, 10, 0.01)
plot(x, pburr(x, alpha=2, rho=-1), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

