\name{QuantMOM}
\alias{QuantMOM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of extreme quantiles using MOM
}
\description{
Compute estimates of an extreme quantile \eqn{Q(1-p)} using the Method of Moments estimates of the EVI.
}
\usage{
QuantMOM(data, gamma, p, plot = FALSE, add = FALSE, 
         main = "Estimates of extreme quantile", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{Moment}}.
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
See Section 4.2.2 of Albrecher et al. (2017) for more details.
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

Dekkers, A.L.M, Einmahl, J.H.J. and de Haan, L. (1989). "A Moment Estimator for the Index of an Extreme-value Distribution." \emph{Annals of Statistics}, 17, 1833--1855.
}
\author{
Tom Reynkens.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ProbMOM}}, \code{\link{Moment}}, \code{\link{QuantGH}}, \code{\link{Quant}}
}
\examples{
data(soa)

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# MOM estimator
M <- Moment(SOAdata)

# Large quantile
p <- 10^(-5)
QuantMOM(SOAdata, p=p, gamma=M$gamma, plot=TRUE)
}
