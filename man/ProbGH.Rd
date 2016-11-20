\name{ProbGH}
\alias{ProbGH}
\alias{ReturnGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and large return periods using generalised Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or large return period \eqn{1/P(X>q)} using the generalised Hill estimates for the EVI.
}
\usage{
ProbGH(data, gamma, q, plot = FALSE, add = FALSE, 
       main = "Estimates of small exceedance probability", ...)

ReturnGH(data, gamma, q, plot = FALSE, add = FALSE, 
         main = "Estimates of large return period", ...)
}
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-2} estimates for the EVI obtained from \code{\link{genHill}}.
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ProbGH} and \code{"Estimates of large return period"} for \code{ReturnGH}.
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
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ProbGH}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{ReturnGH}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Beirlant, J., Vynckier, P. and Teugels, J.L. (1996). "Excess Function and Estimation of the Extreme-value Index". \emph{Bernoulli}, 2, 293--318.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QuantGH}}, \code{\link{genHill}}, \code{\link{ProbMOM}}, \code{\link{Prob}}
}
\examples{
data(soa)

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Hill estimator
H <- Hill(SOAdata)
# Generalised Hill estimator
gH <- genHill(SOAdata, H$gamma)

# Exceedance probability
q <- 10^7
ProbGH(SOAdata, gamma=gH$gamma, q=q, plot=TRUE)

# Return period
q <- 10^7
ReturnGH(SOAdata, gamma=gH$gamma, q=q, plot=TRUE)
}

