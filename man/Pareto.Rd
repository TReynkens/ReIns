\name{Pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Pareto distribution
}
\description{
Density, distribution function, quantile function and random generation for the Pareto distribution (type I).
}
\usage{
dpareto(x, shape, scale = 1, log = FALSE)
ppareto(x, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
qpareto(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
rpareto(n, shape, scale = 1)
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
The shape parameter of the Pareto distribution, a strictly positive number.
}
  \item{scale}{
The scale parameter of the Pareto distribution, a strictly positive number. Its default value is \code{1}.
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
The Cumulative Distribution Function (CDF) of the Pareto distribution is equal to
\eqn{F(x) = 1-(x/scale)^{-shape}} for all \eqn{x \ge scale} and \eqn{F(x)=0} otherwise. Both \code{shape} and \code{scale} need to be strictly positive.
}
\value{
\code{dpareto} gives the density function evaluated in \eqn{x}, \code{ppareto} the CDF evaluated in \code{x} and \code{qpareto} the quantile function evaluated in \eqn{p}. The length of the result is equal to the length of \code{x} or \code{p}.

\code{rpareto} returns a random sample of length \eqn{n}.
}

\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=ptpareto]{tPareto}}, \code{\link[=pgpd]{GPD}},  \code{\link[stats]{Distributions}}
}
\examples{

# Plot of the PDF
x = seq(1,10,0.01)
plot(x, dpareto(x, shape=2), ylab="PDF", type="l")

# Plot of the CDF
x = seq(1,10,0.01)
plot(x, ppareto(x, shape=2), ylab="CDF", type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

