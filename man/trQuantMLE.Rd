\name{trQuantMLE}
\alias{trQuantMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of large quantiles using truncated MLE
}
\description{
This function computes estimates of large quantiles \eqn{Q(1-p)} of the truncated distribution using the ML estimates adapted for upper truncation. Moreover, estimates of large quantiles \eqn{Q_Y(1-p)} of the original distribution \eqn{Y}, which is unobserved, are also computed.
}
\usage{
trQuantMLE(data, gamma, tau, DT, p, Y = FALSE, plot = FALSE, add = FALSE, 
           main = "Estimates of extreme quantile", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
  Vector of \eqn{n} observations. 
}
 \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{trMLE}}.
}
  \item{tau}{
Vector of \eqn{n-1} estimates for the \eqn{\tau} obtained from \code{\link{trMLE}}.
}
\item{DT}{
Vector of \eqn{n-1} estimates for the truncation odds obtained from \code{\link{trDTMLE}}.
}
 \item{p}{
The exceedance probability of the quantile (we estimate \eqn{Q(1-p)} or \eqn{Q_Y(1-p)} for \eqn{p} small).
}
  \item{Y}{
Logical indicating if quantiles from the truncated distribution (\eqn{Q(1-p)}) or from the parent distribution (\eqn{Q_Y(1-p)}) are computed. Default is \code{TRUE}.
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
We observe the truncated r.v. \eqn{X=_d Y | Y<T} where \eqn{T} is the truncation point and \eqn{Y} the untruncated r.v.

Under rough truncation, the quantiles for \eqn{X} are estimated using
\deqn{\hat{Q}_{T,k}(1-p) = X_{n-k,n} +1/(\hat{\tau}_k)([(\hat{D}_{T,k} + (k+1)/(n+1))/(\hat{D}_{T,k}+p)]^{\hat{\xi}_k} -1),}
with \eqn{\hat{\gamma}_k} and \eqn{\hat{\tau}_k} the ML estimates adapted for truncation and \eqn{\hat{D}_T} the estimates for the truncation odds.

The quantiles for \eqn{Y} are estimated using
\deqn{\hat{Q}_{Y,k}(1-p)=X_{n-k,n} +1/(\hat{\tau}_k)([(\hat{D}_{T,k} + (k+1)/(n+1))/(p(\hat{D}_{T,k}+1))]^{\hat{\xi}_k} -1).}

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
Beirlant, J., Fraga Alves, M. I. and Reynkens, T. (2017). "Fitting Tails Affected
by Truncation". \emph{arXiv:1606.02090}.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trMLE}}, \code{\link{trDTMLE}}, \code{\link{trProbMLE}}, \code{\link{trEndpointMLE}}, \code{\link{trTestMLE}}, \code{\link{trQuant}}, \code{\link{Quant}}
}
\examples{
# Sample from GPD truncated at 99\% quantile
gamma <- 0.5
sigma <- 1.5
X <- rtgpd(n=250, gamma=gamma, sigma=sigma, endpoint=qgpd(0.99, gamma=gamma, sigma=sigma))

# Truncated ML estimator
trmle <- trMLE(X, plot=TRUE, ylim=c(0,2))

# Truncation odds
dtmle <- trDTMLE(X, gamma=trmle$gamma, tau=trmle$tau, plot=FALSE)

# Large quantile of X
trQuantMLE(X, gamma=trmle$gamma, tau=trmle$tau, DT=dtmle$DT, plot=TRUE, p=0.005, ylim=c(15,30))

# Large quantile of Y
trQuantMLE(X, gamma=trmle$gamma, tau=trmle$tau, DT=dtmle$DT, plot=TRUE, p=0.005, ylim=c(0,300), 
          Y=TRUE)
}
