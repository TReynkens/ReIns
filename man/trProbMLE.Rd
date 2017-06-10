\name{trProbMLE}
\alias{trProbMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities using truncated MLE
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} using the estimates for the EVI obtained from the ML estimator adapted for upper truncation.
}
\usage{
trProbMLE(data, gamma, tau, DT, q, plot = FALSE, add = FALSE, 
          main = "Estimates of small exceedance probability", ...)
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
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} for \eqn{q} large).
}
 \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The probability is estimated as
\deqn{\hat{p}_{T,k}(q) = (1+ \hat{D}_{T,k}) (k+1)/(n+1) (1+\hat\tau _k(q-X_{n-k,n}))^{-1/\hat{\xi}_k} -\hat{D}_{T,k}}
with \eqn{\hat{\gamma}_k} and  \eqn{\hat{\tau}_k} the ML estimates adapted for truncation and \eqn{\hat{D}_T} the estimates for the truncation odds.

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates.}
\item{q}{The used large quantile.}
}
\references{
Beirlant, J., Fraga Alves, M. I. and Reynkens, T. (2017). "Fitting Tails Affected by Truncation". \emph{Electronic Journal of Statistics}, 11(1), 2026--2065.
}
\author{
Tom Reynkens.
}

\seealso{
\code{\link{trMLE}}, \code{\link{trDTMLE}}, \code{\link{trQuantMLE}}, \code{\link{trEndpointMLE}}, \code{\link{trTestMLE}}, \code{\link{trProb}}, \code{\link{Prob}}
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

# Small exceedance probability
trProbMLE(X, gamma=trmle$gamma, tau=trmle$tau, DT=dtmle$DT, plot=TRUE, q=26, ylim=c(0,0.005))
}