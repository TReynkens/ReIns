\name{ProbEPD}
\alias{ProbEPD}
\alias{ReturnEPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and large return periods using EPD
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or large return period \eqn{1/P(X>q)} using the parameters from the EPD fit.
}
\usage{

ProbEPD(data, q, gamma, kappa, tau, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)

ReturnEPD(data, q, gamma, kappa, tau, plot = FALSE, add = FALSE, 
          main = "Estimates of large return period", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} or \eqn{1/P(X>q)} for \eqn{q} large).
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{EPD}}.
}
  \item{kappa}{
Vector of \eqn{n-1} estimates for \eqn{\kappa} obtained from \code{\link{EPD}}.
}
  \item{tau}{
Vector of \eqn{n-1} estimates for \eqn{\tau} obtained from \code{\link{EPD}}.
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ProbEPD}
and \code{"Estimates of large return period"} for \code{ReturnEPD}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.1 of Albrecher et al. (2017) for more details.
}


\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ProbEPD}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{ReturnEPD}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Joossens, E. and Segers, J. (2009). "Second-Order Refined Peaks-Over-Threshold
Modelling for Heavy-Tailed Distributions." \emph{Journal of Statistical Planning and Inference}, 139, 2800--2815.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{EPD}}, \code{\link{Prob}}
}
\examples{
data(secura)

# EPD estimates for the EVI
epd <- EPD(secura$size, plot=TRUE)

# Compute exceedance probabilities
q <- 10^7
ProbEPD(secura$size, q=q, gamma=epd$gamma, kappa=epd$kappa, tau=epd$tau, plot=TRUE)

# Compute return periods
ReturnEPD(secura$size, q=q, gamma=epd$gamma, kappa=epd$kappa, tau=epd$tau, 
          plot=TRUE, ylim=c(0,10^4))
}
