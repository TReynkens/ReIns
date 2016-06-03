\name{tBurr}
\alias{dtburr}
\alias{ptburr}
\alias{qtburr}
\alias{rtburr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The truncated Burr distribution
}
\description{
Density, distribution function, quantile function and random generation for the truncated Burr distribution (type XII).
}
\usage{
dtburr(x, alpha, rho, eta = 1, endpoint = Inf, log = FALSE)
ptburr(x, alpha, rho, eta = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
qtburr(p, alpha, rho, eta = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
rtburr(n, alpha, rho, eta = 1, endpoint = Inf)
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
The \eqn{\alpha} parameter of the truncated Burr distribution, a strictly positive number.
}
  \item{rho}{
The \eqn{\rho} parameter of the truncated Burr distribution, a strictly negative number.
}
  \item{eta}{
The \eqn{\eta} parameter of the truncated Burr distribution, a strictly positive number.
The default value is \code{1}.
}
  \item{endpoint}{
Endpoint of the truncated Burr distribution. The default value is \code{Inf} for which the truncated Burr distribution corresponds to the ordinary Burr distribution.
}
 \item{log}{
  Logical indicating if the densities are given as \eqn{log(f)}, default is \code{FALSE}.
 }
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as \eqn{log(p)}, default is \code{FALSE}.
 }
}
\details{
The Cumulative Distribution Function (CDF) of the truncated Burr distribution is equal to
\eqn{F_T(x) = F(x) / F(T)} where \eqn{F} is the CDF of the ordinary Burr distribution and \eqn{T} is the endpoint (truncation point) of the truncated Burr distribution.
}
\value{
\code{dtburr} gives the density function evaluated in \eqn{x}, \code{ptburr} the CDF evaluated in \code{x} and \code{qtburr} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rtburr} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=pburr]{Burr}},  \code{\link[stats]{Distributions}}
}
\examples{

# Plot of the PDF
x = seq(0,10,0.01)
plot(x, dtburr(x, alpha=2, rho=-1, endpoint=9), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(0,10,0.01)
plot(x, ptburr(x, alpha=2, rho=-1, endpoint=9), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
