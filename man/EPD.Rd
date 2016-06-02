\name{EPD}
\alias{EPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EPD estimator
}
\description{
Fit the Extended Pareto Distribution (GPD) to the exceedances (peaks) over a threshold. Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
EPD(data, rho = -1, start = NULL, direct = FALSE, warnings = FALSE, 
    logk = FALSE, plot = FALSE, add = FALSE, main = "EPD estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}

  \item{rho}{
A parameter for the \eqn{\rho}-estimator of Fraga Alves et al. (2003)
when strictly positive or choice(s) for \eqn{\rho} if negative. Default is \code{-1}.
}
  \item{start}{
Vector of length 2 containing the starting values for the optimisation. The first element
is the starting value for the estimator of \eqn{\gamma} and the second element is the starting value for the estimator of \eqn{\kappa}. This argument is only used when \code{direct=FALSE}. Default is \code{c(0.1,1)}.
}
  \item{direct}{
Logical indication if the parameters are obtained by directly maximising the log-likelihood function, see Details. Default is \code{TRUE.}
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"EPD estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
We fit the Extended Pareto distribution to the relative excesses over a threshold (X/u).
The EPD has distribution function \eqn{F(x) = 1-(x(1+\kappa-\kappa x^{\tau})^{-1/\gamma}} 
with \eqn{\tau = \rho/\gamma <0<\gamma} and \eqn{\kappa>max(-1,1/\tau)}. 

The parameters are determined using MLE and there are two possible approaches:
maximise the log-likelihood directly (\code{direct=TRUE}) or follow the approach detailed in
Beirlant et al. (2009) (\code{direct=FALSE}). The latter approach uses the score functions of the log-likelihood.

See Section 4.2.1 of Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding estimates for the \eqn{\gamma} parameter of the EPD.}
\item{kappa}{Vector of the corresponding MLE estimates for the \eqn{\kappa} parameter of the EPD.}
\item{tau}{Vector of the corresponding estimates for the \eqn{\tau} parameter of the EPD using Hill estimates and values for \eqn{\rho}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Joossens, E. and Segers, J. (2009). "Second-Order Refined Peaks-Over-Threshold
Modelling for Heavy-Tailed Distributions." \emph{Journal of Statistical Planning and Inference}, 139, 2800--2815.

Fraga Alves, M.I. , Gomes, M.I. and de Haan, L. (2003). "A New Class of Semi-parametric Estimators of the Second Order Parameter." \emph{Portugaliae Mathematica}, 60, 193--214.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GPDmle}}, \code{\link{ReturnEPD}}
}
\examples{
data(secura)

# EPD estimates for the EVI
epd <- EPD(secura$size, plot=TRUE)

# Compute return periods
ReturnEPD(secura$size, 10^10, gamma=epd$gamma, kappa=epd$kappa, 
          tau=epd$tau, plot=TRUE)
}

