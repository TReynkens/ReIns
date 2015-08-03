\name{ES}
\alias{ES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Shortfall
}
\description{
Expected Shortfall (ES) \eqn{ES_p} of the fitted splicing distribution.
}
\usage{
ES(p, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
The probability associated with the ES (we estimate  \eqn{ES_p}).
}
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}} or \code{\link{SpliceFitGPD}}.
}
}
\details{
The Expected Shortfall is defined as
\deqn{ES_p = E(X | X>Q(1-p)) = E(X | X>VaR_p) = VaR_p + 1/p I(VaR_p),}
where \eqn{I(u)=E((X-u)_+)} is the premium of an excess-loss insurance with retention u.

If the CDF is continuous in \eqn{p}, we have \eqn{ES_p=TVaR_p= 1/p \int_0^p VaR_s ds} with
\eqn{TVaR} the Tail Value-at-Risk.

See Section 3.3 of Albrecher et al. (2015) for details.
}
\value{
Vector with the ES associated with each element of p.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang quantiles.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceQuant}}, \code{\link{ExcessSplice}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.6)

p <- seq(0.01, 0.99, 0.01)
# Plot of expected shortfall
plot(p, ES(p, splicefit), type="l", xlab="p", ylab=bquote(ES[p]))
}
}
