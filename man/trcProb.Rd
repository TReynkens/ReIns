\name{trcProb}
\alias{trcProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities using censored truncated Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} using the estimates for the EVI obtained from the Hill estimator adapted for right censoring and upper truncation.
}
\usage{
trcProb(data, censored, gamma1, q, r = 1, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
\item{gamma1}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{trcHill}}.
}
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} for \eqn{q} large).
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
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
\deqn{ (1-km) * ( (q/Z_{n-k,n})^{-1/H_{k,n}^c}-  R_{r,k,n}^{1/\hat{\gamma}_k} ) / (1- R_{r,k,n}^{1/\hat{\gamma}_k})}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{H_{k,n}^c} the Hill estimator adapted for right censoring, \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}, \eqn{R_{r,k,n} = Z_{n-k,n} / Z_{n-r+1,n}}
and  and \eqn{\hat{\gamma}_k} the Hill estimates adapted for right censoring and truncation.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trcHill}}, \code{\link{trciProb}}, \code{\link{trProb}}, \code{\link{cProb}}, \code{\link{KaplanMeier}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=100,shape=shape,endpoint=qpareto(0.99,shape=shape))

# Truncated Hill estimator
trh <- trHill(X,plot=TRUE,ylim=c(0,2))

# Censoring indicator
censored <- sample(0:1,length(X),replace=TRUE)

# Hill estimator adapted for right censoring and 
# upper truncation
trch <- trcHill(X,censored=censored,plot=TRUE,ylim=c(0,2))

# Small probability
trcProb(X,gamma1=trch$gamma1,censored=censored,q=7,plot=FALSE)
}