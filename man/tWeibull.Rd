\name{tWeibull}
\alias{dtweibull}
\alias{ptweibull}
\alias{qtweibull}
\alias{rtweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The truncated Weibull distribution
}
\description{
Density, distribution function, quantile function and random generation for the truncated Weibull distribution.
}
\usage{
dtweibull(x, lambda, tau = 1, endpoint = Inf, log = FALSE)
ptweibull(x, lambda, tau = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
qtweibull(p, lambda, tau = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
rtweibull(n, lambda, tau = 1, endpoint = Inf)
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
  \item{lambda}{
The \eqn{\lambda} parameter of the Weibull distribution, a strictly positive number.
}
  \item{tau}{
The \eqn{\tau} parameter of the Weibull distribution, a strictly positive number, default is 1.
}
  \item{endpoint}{
Endpoint of the truncated Weibull distribution. The default value is \code{Inf} for which the truncated Weibull distribution corresponds to the ordinary Weibull distribution.
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
The Cumulative Distribution Function (CDF) of the truncated Weibull distribution is equal to
\eqn{F_T(x) = F(x) / F(T)} where \eqn{F} is the CDF of the ordinary Weibull distribution and \eqn{T} is the endpoint (truncation point) of the truncated Weibull distribution.
}
\value{
\code{dtweibull} gives the density function evaluated in \eqn{x}, \code{ptweibull} the CDF evaluated in \code{x} and \code{qtweibull} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rtweibull} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Weibull}},  \code{\link[stats]{Distributions}}
}
\examples{

# Plot of the PDF
x = seq(0,10,0.01)
plot(x, dtweibull(x, lambda=2, tau=0.5, endpoint=3), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(0,10,0.01)
plot(x, ptweibull(x, lambda=2, tau=0.5, endpoint=3), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}