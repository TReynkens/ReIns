\name{EPDfit}
\alias{EPDfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit EPD using MLE
}
\description{
Fit the Extended Pareto Distribution (EPD) to data using Maximum Likelihood Estimation (MLE).
}
\usage{
EPDfit(data, tau, start = c(0.1, 1), warnings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{tau}{
 Value for the \eqn{\tau} parameter of the EPD. 
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
See Section 4.2.1 of Albrecher et al. (2017) for more details.
}
\value{
A vector with the MLE estimate for the \eqn{\gamma} parameter of the EPD as the first component and the MLE estimate for the \eqn{\kappa} parameter of the EPD as the second component.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Joossens, E. and Segers, J. (2009). "Second-Order Refined Peaks-Over-Threshold
Modelling for Heavy-Tailed Distributions." \emph{Journal of Statistical Planning and Inference}, 139, 2800--2815.
}
\author{
Tom Reynkens
}


\seealso{
\code{\link{EPD}}, \code{\link{GPDfit}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Fit EPD to last 500 observations
res <- EPDfit(SOAdata/sort(soa$size)[500], tau=-1)
}