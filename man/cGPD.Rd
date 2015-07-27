\name{cGPDmle}
\alias{cGPDmle}
\alias{cPOT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GPD-ML estimator for right censored data
}
\description{
Computes ML estimates of fitting GPD to a peaks over a threshold adapted for right censoring.
}
\usage{
cGPDmle(data, censored, start = c(0.1,1), warnings = FALSE, 
      plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...)
      
cPOT(data, censored, start = c(0.1,1), warnings = FALSE, 
      plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{start}{
Vector of length 2 containing the starting values for the optimisation. The first element
is the starting value for the estimator of \eqn{\gamma_1} and the second element is the starting value for the estimator of \eqn{\sigma_1}. Default is \code{c(0.1,1)}.
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"POT estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The GPD-MLE estimator for the EVI adapted for right censored data is equal to the ordinary GPD-MLE estimator for the EVI divided by the proportion of the \eqn{k} largest observations that is non-censored. The estimates for \eqn{\sigma} are the ordinary GPD-MLE estimates for \eqn{\sigma}.

This estimator is only suitable for right censored data!

\code{cPOT} is the same function but with a different name for compatibility with \code{\link{POT}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma1}{Vector of the corresponding MLE estimates for the \eqn{\gamma_1} parameter of the GPD.}
\item{sigma1}{Vector of the corresponding MLE estimates for the \eqn{\sigma_1} parameter of the GPD.}
}
\references{
Einmahl, J.H.J., Fils-Villetard, A. and Guillou, A. (2008). "Statistics of Extremes Under Random Censoring."
\emph{Bernoulli}, 14, 207--227.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GPDmle}}, \code{\link{cProbGPD}}, \code{\link{cQuantGPD}}, \code{\link{cEPD}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# GPD-ML estimator adapted for right censoring
cpot <- cGPDmle(secura$size,censored=censored,plot=TRUE)
}
