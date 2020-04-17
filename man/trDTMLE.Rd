\name{trDTMLE}
\alias{trDTMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncation odds
}
\description{
Estimates of truncation odds of the truncated probability mass under the untruncated distribution using truncated MLE.
}
\usage{
trDTMLE(data, gamma, tau, plot = FALSE, add = FALSE, main = "Estimates of DT", ...)
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
Logical indicating if the estimates of \eqn{D_T} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{D_T} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of DT"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The truncation odds is defined as
\deqn{D_T=(1-F(T))/F(T)} with \eqn{T} the upper truncation point and \eqn{F} the CDF of the untruncated distribution (e.g. GPD).

We estimate this truncation odds as
\deqn{\hat{D}_T=\max\{ (k+1)/(n+1) ( (1+\hat{\tau}_k E_{1,k})^{-1/\hat{\xi}_k} - 1/(k+1) ) / (1-(1+\hat{\tau}_k E_{1,k})^{-1/\hat{\xi}_k}), 0\}}
with \eqn{E_{1,k} = X_{n,n}-X_{n-k,n}}.

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{DT}{Vector of the corresponding estimates for the truncation odds \eqn{D_T}.}
}
\references{
Beirlant, J., Fraga Alves, M. I. and Reynkens, T. (2017). "Fitting Tails Affected by Truncation". \emph{Electronic Journal of Statistics}, 11(1), 2026--2065.
}
\author{
Tom Reynkens.
}

\seealso{
\code{\link{trMLE}}, \code{\link{trEndpointMLE}}, \code{\link{trProbMLE}}, \code{\link{trQuantMLE}}, \code{\link{trTestMLE}}, \code{\link{trDT}}
}
\examples{
# Sample from GPD truncated at 99\% quantile
gamma <- 0.5
sigma <- 1.5
X <- rtgpd(n=250, gamma=gamma, sigma=sigma, endpoint=qgpd(0.99, gamma=gamma, sigma=sigma))

# Truncated ML estimator
trmle <- trMLE(X, plot=TRUE, ylim=c(0,2))

# Truncation odds
dtmle <- trDTMLE(X, gamma=trmle$gamma, tau=trmle$tau, plot=TRUE, ylim=c(0,0.05))
}