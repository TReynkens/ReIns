\name{QuantGPD}
\alias{QuantGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of extreme quantiles using GPD-MLE
}
\description{
Computes estimates of an extreme quantile \eqn{Q(1-p)} using the GPD fit for the peaks over a threshold.
}
\usage{
QuantGPD(data, gamma, sigma, p, plot = FALSE, add = FALSE, 
         main = "Estimates of extreme quantile", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{GPDmle}}.
}
  \item{sigma}{
Vector of \eqn{n-1} estimates for \eqn{\sigma} obtained from \code{\link{GPDmle}}.
}
  \item{p}{
The exceedance probability of the quantile (we estimate \eqn{Q(1-p)} for \eqn{p} small).
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of extreme quantile"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.2 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens.
}



\seealso{
\code{\link{ProbGPD}}, \code{\link{GPDmle}}, \code{\link{Quant}}
}
\examples{
data(soa)

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# GPD-ML estimator
pot <- GPDmle(SOAdata)

# Large quantile
p <- 10^(-5)
QuantGPD(SOAdata, p=p, gamma=pot$gamma, sigma=pot$sigma, plot=TRUE)
}

