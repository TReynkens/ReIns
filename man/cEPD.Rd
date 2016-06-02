\name{cEPD}
\alias{cEPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EPD estimator for right censored data
}
\description{
Computes the EPD estimates adapted for right censored data.
}
\usage{
cEPD(data, censored, rho = -1, beta = NULL, logk = FALSE, 
     plot = FALSE, add = FALSE, main = "EPD estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{rho}{
A parameter for the \eqn{\rho}-estimator of Fraga Alves et al. (2003)
when strictly positive or choice(s) for \eqn{\rho} if negative. Default is \code{-1}.
}
  \item{beta}{
Parameter for EPD (\eqn{\beta=-\rho/\gamma}). If \code{NULL} (default), \code{beta} is estimated by
\eqn{-\rho/H_{k,n}} with \eqn{H_{k,n}} the Hill estimator.
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
The function \code{EPD} uses \eqn{\tau} which is equal to \eqn{-\beta}.

This estimator is only suitable for right censored data.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma1}{Vector of the corresponding estimates for the \eqn{\gamma} parameter of the EPD.}
\item{kappa1}{Vector of the corresponding MLE estimates for the \eqn{\kappa} parameter of the EPD.}
\item{beta}{Vector values for the \eqn{\beta} parameter of the EPD.}
\item{Delta}{Difference between \code{gamma1} and the Hill estimator for censored data.}
}
\references{
Beirlant, J., Bardoutsos, A., de Wet, T. and Gijbels, I. (2016). "Bias Reduced Tail Estimation for Censored Pareto Type Distributions." \emph{Statistics & Probability Letters}, 109, 78--88.

Fraga Alves, M.I. , Gomes, M.I. and de Haan, L. (2003). "A New Class of Semi-parametric Estimators of the Second Order Parameter." \emph{Portugaliae Mathematica}, 60, 193--214.
}
\author{
Tom Reynkens based on \code{R} code from A. Bardoutsos.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EPD}}, \code{\link{cProbEPD}}, \code{\link{cGPDmle}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# EPD estimator adapted for right censoring
cepd <- cEPD(secura$size,censored=censored,plot=TRUE)
}
