\name{tPareto}
\alias{dtpareto}
\alias{ptpareto}
\alias{qtpareto}
\alias{rtpareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The truncated Pareto distribution
}
\description{
Density, distribution function, quantile function and random generation for the truncated Pareto distribution.
}
\usage{
dtpareto(x, shape, scale = 1, endpoint = Inf, log = FALSE)
ptpareto(x, shape, scale = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
qtpareto(p, shape, scale = 1, endpoint = Inf, lower.tail = TRUE, log.p = FALSE)
rtpareto(n, shape, scale = 1, endpoint = Inf)
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
The shape parameter of the truncated Pareto distribution, a strictly positive number.
}
  \item{scale}{
The scale parameter of the truncated Pareto distribution, a strictly positive number. Its default value is \code{1}.
}
  \item{endpoint}{
Endpoint of the truncated Pareto distribution. The default value is \code{Inf} for which the truncated Pareto distribution corresponds to the ordinary Pareto distribution.
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
The Cumulative Distribution Function (CDF) of the truncated Pareto distribution is equal to
\eqn{F_T(x) = F(x) / F(T)} where \eqn{F} is the CDF of an ordinary Pareto distribution and \eqn{T} is the endpoint (truncation point) of the truncated Pareto distribution.
}
\value{
\code{dtpareto} gives the density function evaluated in \eqn{x}, \code{ptpareto} the CDF evaluated in \code{x} and \code{qtpareto} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rtpareto} returns a random sample of length \eqn{n}.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ppareto]{Pareto}}, \code{\link[stats]{Distributions}}
}
\examples{
# Plot of the PDF
x = seq(1,10,0.01)
plot(x, dtpareto(x, shape=2, endpoint=10), xlab="x", ylab="PDF", type="l")

# Plot of the CDF
x = seq(1,10,0.01)
plot(x, ptpareto(x, shape=2, endpoint=10), xlab="x", ylab="CDF", type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

