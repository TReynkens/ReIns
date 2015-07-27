\name{ciProbEPD}
\alias{ciProbEPD}
\alias{ciReturnEPD}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using interval censored EPD
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the parameters from the EPD fit adapted for interval censoring.
}
\usage{
ciProbEPD(Z, I, censored, gamma1, kappa1, beta, q, plot = FALSE, add = FALSE, 
          main = "Estimates of small exceedance probability", ...)
          
ciReturnEPD(Z, I, censored, gamma1, kappa1, beta, q, plot = FALSE, add = FALSE, 
          main = "Estimates of return period", ...)          
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
  \item{gamma1}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{ciEPD}}.
}
  \item{kappa1}{
Vector of \eqn{n-1} estimates for \eqn{\kappa_1} obtained from \code{\link{ciEPD}}.
}
  \item{beta}{
Vector of \eqn{n-1} estimates for \eqn{\beta} obtained from \code{\link{ciEPD}}.
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ciProbEPD}
and \code{"Estimates of return period"} for \code{ciReturnEPD}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The probability is estimated as \deqn{ (1-tb) * (1-F(x))} with \eqn{F}
the CDF of the EPD with parameters \eqn{\gamma}, \eqn{\kappa} and \eqn{\tau}
and \eqn{tb} the Turnbull estimator for the CDF evaluated in \eqn{Z_{n-k,n}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ciProbEPD}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{cReturnEPD}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ciEPD}}, \code{\link{ProbEPD}}, \code{\link{cProbEPD}}, \code{\link{Prob}}, \code{\link{Turnbull}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Incurred values
I <- secura$size+rchisq(length(secura$size),2)

# EPD estimator adapted for interval censoring
ciepd <- ciEPD(secura$size,I,censored=censored,plot=TRUE)

# Small exceedance proabability
q <- 10^7
ciProbEPD(secura$size,I,censored=censored,gamma1=ciepd$gamma1,
          kappa1=ciepd$kappa1,beta=ciepd$beta,q=q,plot=TRUE)   
          
# Return period
ciReturnEPD(secura$size,I,censored=censored,gamma1=ciepd$gamma1,
            kappa1=ciepd$kappa1,beta=ciepd$beta,q=q,plot=TRUE)          
}