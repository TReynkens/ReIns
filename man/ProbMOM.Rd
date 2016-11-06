\name{ProbMOM}
\alias{ProbMOM}
\alias{ReturnMOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using MOM
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the Method of Moments estimates for the EVI.
}
\usage{
ProbMOM(data, gamma, q, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)

ReturnMOM(data, gamma, q, plot = FALSE, add = FALSE, 
          main = "Estimates of return period", ...)
}
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{Moment}}.
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ProbMOM}
and \code{"Estimates of return period"} for \code{ReturnMOM}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.2 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ProbMOM}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{ReturnMOM}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Dekkers, A.L.M, Einmahl, J.H.J. and de Haan, L. (1989). "A Moment Estimator for the Index of an Extreme-value Distribution." \emph{Annals of Statistics}, 17, 1833--1855.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QuantMOM}}, \code{\link{Moment}}, \code{\link{ProbGH}}, \code{\link{Prob}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# MOM estimator
M <- Moment(SOAdata)

# Exceedance probability
q = 10^7
ProbMOM(SOAdata, gamma=M$gamma, q=q, plot=TRUE)

# Return period
q = 10^7
ReturnMOM(SOAdata, gamma=M$gamma, q=q, plot=TRUE)
}

