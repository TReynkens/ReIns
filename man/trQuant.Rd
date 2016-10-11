\name{trQuant}
\alias{trQuant}
\alias{trQuantW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of large quantiles using truncated Hill
}
\description{
\code{trQuant} computes estimates of large quantiles \eqn{Q(1-p)} of the truncated data using the estimates for the EVI obtained from the Hill estimator adapted for upper truncation. \code{trQuantW} computes estimates of large quantiles \eqn{Q_W(1-p)} of the original data \eqn{W} which is unobserved.
}
\usage{
trQuant(data, r = 1, rough = TRUE, gamma, DT, p, plot = FALSE, add = FALSE, 
        main = "Estimates of extreme quantile", ...)
        
trQuantW(data, gamma, DT, p, plot = FALSE, add = FALSE, 
        main = "Estimates of extreme quantile", ...)        
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations (truncated data). 
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
}
  \item{rough}{
Logical indicating if rough truncation is present, default is \code{TRUE}.
}
\item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{trHill}}.
}
\item{DT}{
Vector of \eqn{n-1} estimates for the odds ratio obtained from \code{\link{trDT}}.
}
  \item{p}{
The exceedance probability of the quantile (we estimate \eqn{Q(1-p)} for \eqn{p} small).
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
We observe the truncated r.v. \eqn{X=W | W<T} where \eqn{T} is the truncation point and \eqn{W} the untruncated r.v..

Under rough truncation, the quantiles for \eqn{X} are estimated using
\deqn{\hat{Q}(1-p)=X_{n-k,n}  ((\hat{D}_T + (k+1)/(n+1))/(\hat{D}_T+p))^{\hat{\gamma}_k},}
with \eqn{\hat{\gamma}_k} the Hill estimates adapted for truncation and \eqn{\hat{D}_T} the estimates for the odds ratio.

Under light truncation, the quantiles are estimated using the Weissman estimator with the Hill estimates replaced
by the truncated Hill estimates:
\deqn{\hat{Q}(1-p)=X_{n-k,n}  ((k+1)/((n+1)p))^{\hat{\gamma}_k}.}

To decide between light and rough truncation, one can use the test implemented in \code{\link{trTest}}.

The quantiles for \eqn{W} are estimated using
\deqn{\hat{Q}_W(1-p)=X_{n-k,n} ( (\hat{D}_T + (k+1)/(n+1)) / (p(1+\hat{D}_T))^{\hat{\gamma}_k}.}

See Beirlant et al. (2016) or Section 4.2.3 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Fraga Alves, M.I. and Gomes, M.I. (2016). "Tail fitting for Truncated
Pareto-type Distributions." \emph{Extremes}, accepted.
}
\author{
Tom Reynkens based on \code{R} code of Dries Cornilly.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trHill}}, \code{\link{trDT}}, \code{\link{trProb}}, \code{\link{trEndpoint}}, \code{\link{trTest}}, \code{\link{Quant}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000, shape=shape, endpoint=qpareto(0.99, shape=shape))

# Truncated Hill estimator
trh <- trHill(X, plot=TRUE, ylim=c(0,2))

# Odds ratio
dt <- trDT(X, gamma=trh$gamma, plot=TRUE, ylim=c(0,2))

# Large quantile
p <- 10^(-5)
# Truncated distribution
trQuant(X, gamma=trh$gamma, DT=dt$DT, p=p, plot=TRUE)
# Original distribution
trQuantW(X, gamma=trh$gamma, DT=dt$DT, p=p, plot=TRUE, ylim=c(0,1000))

}