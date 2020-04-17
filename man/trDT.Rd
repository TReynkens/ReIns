\name{trDT}
\alias{trDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncation odds
}
\description{
Estimates of truncation odds of the truncated probability mass under the untruncated distribution using truncated Hill.
}
\usage{
trDT(data, r = 1, gamma, plot = FALSE, add = FALSE, main = "Estimates of DT", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
}
\item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{trHill}}.
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
\deqn{D_T=(1-F(T))/F(T)} with \eqn{T} the upper truncation point and \eqn{F} the CDF of the untruncated distribution (e.g. Pareto distribution).

We estimate this truncation odds as
\deqn{\hat{D}_T=\max\{ (k+1)/(n+1) ( R_{r,k,n}^{1/\gamma_k} - 1/(k+1) ) / (1-R_{r,k,n}^{1/\gamma_k}), 0\}}
with \eqn{R_{r,k,n} = X_{n-k,n} / X_{n-r+1,n}}.

See Beirlant et al. (2016) or Section 4.2.3 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{DT}{Vector of the corresponding estimates for the truncation odds \eqn{D_T}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Fraga Alves, M.I. and Gomes, M.I. (2016). "Tail fitting for Truncated and Non-truncated Pareto-type Distributions." \emph{Extremes}, 19, 429--462.
}
\author{
Tom Reynkens based on \code{R} code of Dries Cornilly.
}

\seealso{
\code{\link{trHill}}, \code{\link{trEndpoint}}, \code{\link{trQuant}}, \code{\link{trDTMLE}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000, shape=shape, endpoint=qpareto(0.99, shape=shape))

# Truncated Hill estimator
trh <- trHill(X, plot=TRUE, ylim=c(0,2))

# Truncation odds
dt <- trDT(X, gamma=trh$gamma, plot=TRUE, ylim=c(0,0.05))
}