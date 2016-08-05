\name{ES}
\alias{ES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Shortfall
}
\description{
Compute Expected Shortfall (ES) \eqn{ES_p} of the fitted spliced distribution.
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
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
}
\details{
The Expected Shortfall is defined as
\deqn{ES_p = E(X | X>Q(1-p)) = E(X | X>VaR_p) = VaR_p + 1/p I(VaR_p),}
where \eqn{I(u)=E((X-u)_+)} is the premium of an excess-loss insurance with retention u.

If the CDF is continuous in \eqn{p}, we have \eqn{ES_p=TVaR_p= 1/p \int_0^p VaR_s ds} with
\eqn{TVaR} the Tail Value-at-Risk.

See Reynkens et al. (2016) and Section 4.5 of Albrecher et al. (2016) for more details.
}
\value{
Vector with the ES corresponding to each element of \eqn{p}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2016). "Modeling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang quantiles.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qSplice}}, \code{\link{ExcessSplice}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}},  \code{\link{SpliceFitGPD}} 
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)

p <- seq(0.01, 0.99, 0.01)
# Plot of expected shortfall
plot(p, ES(p, splicefit), type="l", xlab="p", ylab=bquote(ES[p]))
}
}
