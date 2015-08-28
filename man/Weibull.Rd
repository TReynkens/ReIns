\name{Weibull}
\alias{dweibull}
\alias{pweibull}
\alias{qweibull}
\alias{rweibull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Weibull distribution
}
\description{
Density, distribution function, quantile function and random generation for the Weibull distribution.
}
\usage{
dweibull(x, lambda, tau = 1, log = FALSE)
pweibull(x, lambda, tau = 1, lower.tail = TRUE, log.p = FALSE)
qweibull(p, lambda, tau = 1, lower.tail = TRUE, log.p = FALSE)
rweibull(n, lambda, tau = 1)
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
The Cumulative Distribution Function (CDF) of the Weibull distribution is equal to
\eqn{F(x) = 1-exp(-\lambda*x^{\tau})} for all \eqn{x \ge 0} and \eqn{F(x)=0} otherwise. We need that \eqn{\lambda>0} and \eqn{\tau>0}. 
}
\value{
\code{dweibull} gives the density function evaluated in \eqn{x}, \code{pweibull} the CDF evaluated in \code{x} and \code{qweibull} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rweibull} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{Exponential}}, \code{\link[=dtweibull]{tWeibull}},  \code{\link[stats]{Distributions}}
}
\examples{

# Plot of the PDF
x = seq(0,10,0.01)
plot(x, dweibull(x, lambda=2, tau=0.5), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(0,10,0.01)
plot(x, pweibull(x, lambda=2, tau=0.5), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}