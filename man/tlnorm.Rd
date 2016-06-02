\name{tlnorm}
\alias{dtlnorm}
\alias{ptlnorm}
\alias{qtlnorm}
\alias{rtlnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The truncated log-normal distribution
}
\description{
Density, distribution function, quantile function and random generation for the truncated log-normal distribution.
}
\usage{
dtlnorm(x, meanlog = 0, sdlog = 1, endpoint = Inf, log = FALSE)
ptlnorm(x, meanlog = 0, sdlog = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
qtlnorm(p, meanlog = 0, sdlog = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
rtlnorm(n, meanlog = 0, sdlog = 1, endpoint = Inf)
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
  \item{meanlog}{
Mean of the distribution on the log scale, default is 0.
}
  \item{sdlog}{
Standard deviation of the distribution on the log scale, default is 1.
}
  \item{endpoint}{
Endpoint of the truncated log-normal distribution. The default value is \code{Inf} for which the truncated log-normal distribution corresponds to the ordinary log-normal distribution.
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
The Cumulative Distribution Function (CDF) of the truncated log-normal distribution is equal to
\eqn{F_T(x) = F(x) / F(T)} where \eqn{F} is the CDF of the ordinary log-normal distribution and \eqn{T} is the endpoint (truncation point) of the truncated log-normal distribution.
}
\value{
\code{dtlnorm} gives the density function evaluated in \eqn{x}, \code{ptlnorm} the CDF evaluated in \code{x} and \code{qtlnorm} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rtlnorm} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{Lognormal}},  \code{\link[stats]{Distributions}}
}
\examples{

# Plot of the PDF
x = seq(0,10,0.01)
plot(x, dtlnorm(x, endpoint=9), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(0,10,0.01)
plot(x, ptlnorm(x, endpoint=9), xlab="x", ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
