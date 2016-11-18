\name{GPDfit}
\alias{GPDfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit GPD using MLE
}
\description{
Fit the Generalised Pareto Distribution (GPD) to data using Maximum Likelihood Estimation (MLE).
}
\usage{
GPDfit(data, start = c(0.1, 1), warnings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{start}{
Vector of length 2 containing the starting values for the optimisation. The first element
is the starting value for the estimator of \eqn{\gamma} and the second element is the starting value for the estimator of \eqn{\sigma}. Default is \code{c(0.1,1)}.
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
}
\details{
See Section 4.2.2 in Albrecher et al. (2017) for more details.
}
\value{
A vector with the MLE estimate for the \eqn{\gamma} parameter of the GPD as the first component and the MLE estimate for the \eqn{\sigma} parameter of the GPD as the second component.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur and \code{R} code from Klaus Herrmann.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GPDmle}}, \code{\link{EPDfit}}
}
\examples{
data(soa)

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Fit GPD to last 500 observations
res <- GPDfit(SOAdata-sort(soa$size)[500])
}
