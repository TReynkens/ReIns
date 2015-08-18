\name{qSplice}
\alias{qSplice}
\alias{VaR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantiles of splicing fit
}
\description{
Compute quantiles (Value-at-Risk) of the fitted splicing distribution.
}
\usage{
qSplice(p, splicefit, lower.tail = TRUE, log.p = FALSE)

VaR(p, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
The probability of the quantile for \code{SpliceQuant} (we estimate \eqn{Q(p)}) and the exceedance probability for \code{VaR} (we estimate \eqn{VaR_p=Q(1-p)}).
}
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}} or \code{\link{SpliceFitGPD}}.
}
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as log(p), default is \code{FALSE}.
 }
}
\details{ 
See Section 3.3 of Albrecher et al. (2015) for details.

Note that \code{VaR(p, splicefit)} corresponds to \code{qSplice(p, splicefit, log = FALSE, lower.tail = FALSE)}.
}
\value{
Vector of quantiles \eqn{Q(p)} for \code{SpliceQuant} and vector of quantiles \eqn{VaR_p=Q(1-p)} for \code{VaR}.
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
\code{\link{rSplice}}, \code{\link{pSplice}}, \code{\link{dSplice}}, \code{\link{ES}}, \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}}
}

\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.6)

p <- seq(0,1,0.01)

# Plot of quantiles
plot(p, qSplice(p, splicefit), type="l", xlab="p", ylab="Q(p)")

# Plot of VaR
plot(p, VaR(p, splicefit), type="l", xlab="p", ylab=bquote(VaR[p]))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }

