\name{CTE}
\alias{CTE}
\alias{ES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional Tail Expectation
}
\description{
Compute Conditional Tail Expectation (CTE) \eqn{CTE_{1-p}} of the fitted spliced distribution.
}
\usage{
CTE(p, splicefit)

ES(p, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
The probability associated with the CTE (we estimate  \eqn{CTE_{1-p}}).
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
}
\details{
The Conditional Tail Expectation is defined as
\deqn{CTE_{1-p} = E(X | X>Q(1-p)) = E(X | X>VaR_{1-p}) = VaR_{1-p} + \Pi(VaR_{1-p})/p,}
where \eqn{\Pi(u)=E((X-u)_+)} is the premium of the excess-loss insurance with retention u.

If the CDF is continuous in \eqn{p}, we have \eqn{CTE_{1-p}=TVaR_{1-p}= 1/p \int_0^p VaR_{1-s} ds} with
\eqn{TVaR} the Tail Value-at-Risk.

See Reynkens et al. (2017) and Section 4.6 of Albrecher et al. (2017) for more details.

The \code{ES} function is the same function as \code{CTE} but is deprecated. 
}
\value{
Vector with the CTE corresponding to each element of \eqn{p}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{Insurance: Mathematics and Economics}, accepted. Available on arXiv:1608.01566.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the mixed Erlang quantiles.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qSplice}}, \code{\link{ExcessSplice}}, \code{\link{SpliceFit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}},  \code{\link{SpliceFitGPD}} 
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)

p <- seq(0.01, 0.99, 0.01)
# Plot of CTE
plot(p, CTE(p, splicefit), type="l", xlab="p", ylab=bquote(CTE[1-p]))
}
}
