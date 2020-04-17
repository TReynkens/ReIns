\name{genHill}
\alias{genHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalised Hill estimator
}
\description{
Computes the generalised Hill estimator for real extreme value indices as a function of the tail parameter \eqn{k}.
Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
genHill(data, gamma, logk = FALSE, plot = FALSE, add = FALSE, 
        main = "Generalised Hill estimates of the EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, typically Hill estimates are used.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{\log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Generalised Hill estimates of the EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}

\details{
The generalised Hill estimator is an estimator for the slope of the \eqn{k} last points of the generalised QQ-plot:
\deqn{\hat{\gamma}^{GH}_{k,n}=1/k\sum_{j=1}^k \log UH_{j,n}- \log UH_{k+1,n}}
with \eqn{UH_{j,n}=X_{n-j,n}H_{j,n}} the UH scores and \eqn{H_{j,n}} the Hill estimates.
This is analogous to the (ordinary) Hill estimator which is the estimator of the slope of the \eqn{k} last points of the Pareto QQ-plot when using constrained least squares. 

See Section 4.2.2 of Albrecher et al. (2017) for more details.
}

\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding generalised Hill estimates.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Beirlant, J., Vynckier, P. and Teugels, J.L. (1996). "Excess Function and Estimation of the Extreme-value Index". \emph{Bernoulli}, 2, 293--318.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{genQQ}}, \code{\link{Moment}}
}
\examples{
data(soa)

# Hill estimator
H <- Hill(soa$size, plot=FALSE)
# Moment estimator
M <- Moment(soa$size)
# Generalised Hill estimator
gH <- genHill(soa$size, gamma=H$gamma)

# Plot estimates
plot(H$k[1:5000], M$gamma[1:5000], xlab="k", ylab=expression(gamma), type="l", ylim=c(0.2,0.5))
lines(H$k[1:5000], gH$gamma[1:5000], lty=2)
legend("topright", c("Moment", "Generalised Hill"), lty=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

