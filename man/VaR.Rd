\name{VaR}
\alias{VaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
VaR of splicing fit
}
\description{
Compute Value-at-Risk (\eqn{VaR_{1-p}=Q(1-p)}) of the fitted spliced distribution.
}
\usage{
VaR(p, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
The exceedance probability (we estimate \eqn{VaR_{1-p}=Q(1-p)}).
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
}
\details{ 
See Reynkens et al. (2016) and Section 4.6 of Albrecher et al. (2017) for details.

Note that \code{VaR(p, splicefit)} corresponds to \code{qSplice(p, splicefit, lower.tail = FALSE)}.
}
\value{
Vector of quantiles \eqn{VaR_{1-p}=Q(1-p)}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

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
\code{\link{qSplice}}, \code{\link{ES}}, \code{\link{SpliceFit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceFitGPD}}
}

\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)



p <- seq(0,1,0.01)

# Plot of quantiles
plot(p, qSplice(p, splicefit), type="l", xlab="p", ylab="Q(p)")

# Plot of VaR
plot(p, VaR(p, splicefit), type="l", xlab="p", ylab=bquote(VaR[1-p]))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


