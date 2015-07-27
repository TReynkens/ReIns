\name{ciEPD}
\alias{ciEPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EPD estimator for interval censored data
}
\description{
EPD estimator adapted for interval censored data.
}
\usage{
ciEPD(Z, I, censored, threshold = NULL, rho = -1, beta = NULL, 
      start = NULL, warnings = FALSE, plot = FALSE, add = FALSE, 
      main = "EPD estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals. In an 
insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{threshold}{
A vector of thresholds. If \code{NULL} (default), the order statistics \eqn{Z_{n-k,n}}
for \eqn{k=1,...,n-1} are used.
}
  \item{rho}{
A parameter for the \eqn{\rho}-estimator of Fraga Alves et al. (2003)
when strictly positive or choice(s) for \eqn{\rho} if negative. Default is \code{-1}.
}
  \item{beta}{
Parameter for EPD (\eqn{\beta=-\tau/\gamma}). If \code{NULL} (default), \code{beta} is estimated by
\eqn{-\rho/H_{k,n}} with \eqn{H_{k,n}} the Hill estimator.
}
  \item{start}{
Vector of length 2 containing the starting values for the optimisation. The first element
is the starting value for the estimator of \eqn{\gamma_1} and the second element is the starting value for the estimator of \eqn{\kappa_1}. When \code{NULL} (default), we use the Hill estimator (\code{threshold=NULL}) or 1 (otherwise) as starting value for \eqn{\gamma_1} and 0 for \eqn{\kappa_1}.
}
 \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates of \eqn{\gamma_1} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma_1} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"EPD estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\value{
A list with following components:
\item{gamma1}{Vector of the corresponding estimates for the \eqn{\gamma} parameter of the EPD.}
\item{kappa1}{Vector of the corresponding MLE estimates for the \eqn{\kappa} parameter of the EPD.}
\item{beta}{Vector values for the \eqn{\beta} parameter of the EPD.}
\item{k}{Vector of the values of the tail parameter \eqn{k}, only returned when \code{threshold=NULL}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EPD}}, \code{\link{cEPD}}, \code{\link{ciProbEPD}}, \code{\link{ciHill}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Incurred values
I <- secura$size+rchisq(length(secura$size),2)

# EPD estimator adapted for interval censoring
ciepd <- ciEPD(secura$size,I,censored=censored,plot=TRUE)
}
