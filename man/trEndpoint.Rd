\name{trEndpoint}
\alias{trEndpoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of endpoint
}
\description{
Estimator of endpoint using truncated Hill estimates.
}
\usage{
trEndpoint(data, r = 1, gamma, DT, plot = FALSE, add = FALSE,
          main = "Estimates of Endpoint", ...)
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
\item{DT}{
Vector of \eqn{n-1} estimates for the odds ratio obtained from \code{\link{trDT}}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of DT"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The endpoint is estimated as
\deqn{\log(\hat{T}_{k,n}) = \max( \log(X_{n-k,n}) + \hat{\gamma}_k \log(1+(k+1)/((n+1)\hat{D_T})), \log(X_{n,n}) ),}
with \eqn{\hat{\gamma}_k} the Hill estimates adapted for truncation and \eqn{hat{D_T}} the estimates for the odds ratio.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Tk}{Vector of the corresponding estimates for the endpoint \eqn{T}.}
}
\references{
Beirlant, J., Fraga Alves, M.I., Gomes, M.I. and Meerschaert, M.M. (2015). "Extreme Value Statistics for Truncated
Pareto-type Distributions", \emph{http://arxiv.org/abs/1410.4097}.
}
\author{
Tom Reynkens based on \code{R} code of Dries Cornilly.
}

\seealso{
\code{\link{trHill}}, \code{\link{trDT}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000,shape=shape,endpoint=qpareto(0.99,shape=shape))

# Truncated Hill estimator
trh <- trHill(X,plot=TRUE,ylim=c(0,2))

# Odds ratio
dt <- trDT(X,gamma=trh$gamma,plot=TRUE,ylim=c(0,2))

# Endpoint
trEndpoint(X,gamma=trh$gamma, DT=dt$DT,plot=TRUE,ylim=c(8,12))
abline(h=qpareto(0.99,shape=shape),lty=2)
}
