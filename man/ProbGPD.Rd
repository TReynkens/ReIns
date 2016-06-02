\name{ProbGPD}
\alias{ProbGPD}
\alias{ReturnGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using GPD-MLE
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the GPD fit for the peaks over a threshold.
}
\usage{
ProbGPD(data, gamma, sigma, q, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)

ReturnGPD(data, gamma, sigma, q, plot = FALSE, add = FALSE, 
        main = "Estimates of return period", ...)
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
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} or \eqn{1/P(X>q)} for \eqn{q} large).
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ProbGPD}
and \code{"Estimates of return period"} for \code{ReturnGPD}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.2 in Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ProbGPD}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{ReturnGPD}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QuantGPD}}, \code{\link{GPDmle}}, \code{\link{Prob}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# GPD-ML estimator
pot <- GPDmle(SOAdata)

# Exceedance probability
q = 10^7
ProbGPD(SOAdata, gamma=pot$gamma, sigma=pot$sigma, q=q, plot=TRUE)

# Return period
q = 10^7
ReturnGPD(SOAdata, gamma=pot$gamma, sigma=pot$sigma, q=q, plot=TRUE)
}

