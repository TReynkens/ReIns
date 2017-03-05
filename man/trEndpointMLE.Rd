\name{trEndpointMLE}
\alias{trEndpointMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of endpoint
}
\description{
Estimator of endpoint using truncated ML estimates.
}
\usage{
trEndpointMLE(data, gamma, tau, plot = FALSE, add = FALSE, 
              main = "Estimates of endpoint", ...)
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
   \item{plot}{
Logical indicating if the estimates of \eqn{T} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{T} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of endpoint"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The endpoint is estimated as
\deqn{\hat{T}_{k} = X_{n-k,n} + 1/\hat{\tau}_k[( (1-1/k)/((1+ \hat{\tau}_k (X_{n,n}-X_{n-k,n}))^{-1/\hat{\xi}_k}-1/k))^{\hat{\xi}_k} -1]}
with \eqn{\hat{\gamma}_k} and \eqn{\hat{\tau}_k} the truncated ML estimates for \eqn{\gamma} and \eqn{\tau}.

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Tk}{Vector of the corresponding estimates for the endpoint \eqn{T}.}
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
\code{\link{trMLE}}, \code{\link{trDTMLE}}, \code{\link{trProbMLE}}, \code{\link{trQuantMLE}}, \code{\link{trTestMLE}}, \code{\link{trEndpoint}}
}
\examples{
# Sample from GPD truncated at 99\% quantile
gamma <- 0.5
sigma <- 1.5
X <- rtgpd(n=250, gamma=gamma, sigma=sigma, endpoint=qgpd(0.99, gamma=gamma, sigma=sigma))

# Truncated ML estimator
trmle <- trMLE(X, plot=TRUE, ylim=c(0,2))

# Endpoint
trEndpointMLE(X, gamma=trmle$gamma, tau=trmle$tau, plot=TRUE, ylim=c(0,50))
}
