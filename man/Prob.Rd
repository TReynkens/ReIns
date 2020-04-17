\name{Prob}
\alias{Prob}
\alias{Return}
\alias{Weissman.p}
\alias{Weissman.r}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weissman estimator of small exceedance probabilities and large return periods
}
\description{
Compute estimates of a small exceedance probability \eqn{P(X>q)} or large return period \eqn{1/P(X>q)} using the approach of Weissman (1978).
}
\usage{
Prob(data, gamma, q, plot = FALSE, add = FALSE, 
     main = "Estimates of small exceedance probability", ...)
          
Return(data, gamma, q, plot = FALSE, add = FALSE, 
       main = "Estimates of large return period", ...)
        
Weissman.p(data, gamma, q, plot = FALSE, add = FALSE, 
           main = "Estimates of small exceedance probability", ...)
          
Weissman.r(data, gamma, q, plot = FALSE, add = FALSE, 
           main = "Estimates of large return period", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, typically Hill estimates are used.
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
Title for the plot, default is \code{"Estimates of extreme quantile"} for \code{Prob}
and \code{"Estimates of large return period"} for \code{Return}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
See Section 4.2.1 of Albrecher et al. (2017) for more details.

\code{Weissman.p} and \code{Weissman.r} are the same functions as \code{Prob} and \code{Return} but with a different name for compatibility with the old \code{S-Plus} code. 
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{Prob}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{Return}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Weissman, I. (1978). "Estimation of Parameters and Large Quantiles Based on the \emph{k} Largest Observations." \emph{Journal of the American Statistical Association}, 73, 812--815.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Quant}}
}
\examples{
data(soa)

# Look at last 500 observations of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Hill estimator
H <- Hill(SOAdata)

# Exceedance probability
q <- 10^6
# Weissman estimator
Prob(SOAdata,gamma=H$gamma,q=q,plot=TRUE)

# Return period
q <- 10^6
# Weissman estimator
Return(SOAdata,gamma=H$gamma,q=q,plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

