\name{trTestMLE}
\alias{trTestMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for truncated GPD tails
}
\description{
Test between non-truncated GPD tails \emph{(light truncation)} and truncated GPD tails \emph{(rough truncation)}.
}
\usage{
trTestMLE(data, gamma, tau, alpha = 0.05, plot = TRUE, main = "Test for truncation", ...)
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
 \item{alpha}{
 The used significance level, default is \code{0.05}.
}
 \item{plot}{
Logical indicating if the P-values should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Test for truncation"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
We want to test
\eqn{H_0: X} has non-truncated GPD tails vs.
\eqn{H_1: X} has truncated GPD tails.
Let \eqn{\hat{\gamma}_k} and \eqn{\hat{\tau}_k} be the truncated MLE estimates for \eqn{\gamma} and \eqn{\tau}.
The test statistic is then
\deqn{T_{k,n}=k (1+\hat{\tau} (X_{n,n}-X_{-k,n}) )^{-1/\hat{\xi}_k}}
which is asymptotically standard exponentially distributed.
We reject \eqn{H_0} on level \eqn{\alpha} if 
\eqn{T_{k,n}>\ln (1/{\alpha)}}. The corresponding P-value is given by
\eqn{\exp(-T_{k,n})}.

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{testVal}{Corresponding test values.}
\item{critVal}{Critical value used for the test, i.e. \eqn{\ln(1/\alpha)}.}
\item{Pval}{Corresponding P-values.}
\item{Reject}{Logical vector indicating if the null hypothesis is rejected for a certain value of \code{k}.}
}
\references{
Beirlant, J., Fraga Alves, M. I. and Reynkens, T. (2017). "Fitting Tails Affected by Truncation". \emph{Electronic Journal of Statistics}, 11(1), 2026--2065.
}
\author{
Tom Reynkens.
}

\seealso{
\code{\link{trMLE}}, \code{\link{trDTMLE}}, \code{\link{trProbMLE}}, \code{\link{trEndpointMLE}}, \code{\link{trTestMLE}}, \code{\link{trTest}}
}
\examples{
# Sample from GPD truncated at 99\% quantile
gamma <- 0.5
sigma <- 1.5
X <- rtgpd(n=250, gamma=gamma, sigma=sigma, endpoint=qgpd(0.99, gamma=gamma, sigma=sigma))

# Truncated ML estimator
trmle <- trMLE(X, plot=TRUE, ylim=c(0,2))

# Test for truncation
trTestMLE(X, gamma=trmle$gamma, tau=trmle$tau)
}