\name{QuantGH}
\alias{QuantGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of extreme quantiles using generalised Hill
}
\description{
Compute estimates of an extreme quantile \eqn{Q(1-p)} using generalised Hill estimates of the EVI.
}
\usage{
QuantGH(data, gamma, p, plot = FALSE, add = FALSE, 
        main = "Estimates of extreme quantile", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-2} estimates for the EVI obtained from \code{\link{genHill}}.
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
See Section 4.2.2 of Albrecher et al. (2016) for more details.
}

\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Beirlant, J., Vynckier, P. and Teugels, J.L. (1996). "Excess Function and Estimation of the Extreme-value Index". \emph{Bernoulli}, 2, 293--318.
}
\author{
Tom Reynkens.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ProbGH}}, \code{\link{genHill}}, \code{\link{QuantMOM}}, \code{\link{Quant}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Hill estimator
H <- Hill(SOAdata)
# Generalised Hill estimator
gH <- genHill(SOAdata, H$gamma)

# Large quantile
p <- 10^(-5)
QuantGH(SOAdata, p=p, gamma=gH$gamma, plot=TRUE)
}
