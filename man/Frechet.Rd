\name{Fréchet}
\alias{dfrechet}
\alias{pfrechet}
\alias{qfrechet}
\alias{rfrechet}
\encoding{UTF-8}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Frechet distribution
}
\description{
Density, distribution function, quantile function and random generation for the Fréchet distribution (inverse Weibull distribution).
}
\usage{
dfrechet(x, shape, loc = 0, scale = 1, log = FALSE)
pfrechet(x, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qfrechet(p, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rfrechet(n, shape, loc = 0, scale = 1)
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
  \item{shape}{
Shape parameter of the Fréchet distribution.
}
  \item{loc}{
Location parameter of the Fréchet distribution, default is 0.
}
  \item{scale}{
Scale parameter of the Fréchet distribution, default is 1.
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
The Cumulative Distribution Function (CDF) of the Fréchet distribution is equal to
\eqn{F(x) = exp(-((x-loc)/scale)^{-shape})} for all \eqn{x \ge loc} and \eqn{F(x)=0} otherwise. Both \code{shape} and \code{scale} need to be strictly positive.
}
\value{
\code{dfrechet} gives the density function evaluated in \eqn{x}, \code{pfrechet} the CDF evaluated in \code{x} and \code{qfrechet} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rfrechet} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ptfrechet]{tFréchet}},  \code{\link[stats]{Distributions}}
}
\examples{
# Plot of the PDF
x <- seq(1,10,0.01)
plot(x, dfrechet(x, shape=2), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x <- seq(1,10,0.01)
plot(x, pfrechet(x, shape=2), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

