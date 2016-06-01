\name{Quant}
\alias{Quant}
\alias{Weissman.q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weissman estimator of extreme quantiles
}
\description{
Compute estimates of an extreme quantile \eqn{Q(1-p)} using the approach of Weissman (1978).
}
\usage{
Quant(data, gamma, p, plot = FALSE, add = FALSE, 
            main = "Estimates of extreme quantile", ...)
            
Weissman.q(data, gamma, p, plot = FALSE, add = FALSE, 
            main = "Estimates of extreme quantile", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, typically Hill estimates are used.
}
  \item{p}{
The exceedance probability of the quantile (we estimate \eqn{Q(1-p)} for \eqn{p} small).
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates as a function of \eqn{k} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of extreme quantile"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.1 of Albrecher et al. (2016) for more details.

\code{Weissman.q} is the same function but with a different name for compatibility with the old \code{S-Plus} code. 
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley.

Weissman, I. (1978). "Estimation of Parameters and Large Quantiles Based on the \emph{k} Largest Observations." \emph{Journal of the American Statistical Association}, 73, 812--815.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Prob}}, \code{\link{Quant.2oQV}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata = sort(soa$size)[length(soa$size)-(0:499)]

# Hill estimator
H <- Hill(SOAdata)
# Bias-reduced estimator (QV)
H_QV <- Hill.2oQV(SOAdata)

# Exceedance probability
p <- 10^(-5)
# Weissman estimator
Quant(SOAdata,gamma=H$gamma,p=p,plot=TRUE)

# Second order Weissman estimator (QV)
Quant.2oQV(SOAdata,gamma=H_QV$gamma,beta=H_QV$beta,b=H_QV$b,p=p,
                add=TRUE,lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
