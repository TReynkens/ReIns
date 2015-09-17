\name{trProb}
\alias{trProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities using truncated Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} using the estimates for the EVI obtained from the Hill estimator adapted for upper truncation.
}
\usage{
trProb(data, r = 1, gamma, q, warnings = TRUE, plot = FALSE, add = FALSE, 
       main = "Estimates of small exceedance probability", ...)
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
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} for \eqn{q} large).
}
  \item{warnings}{
Logical indicating if warnings are shown, default is \code{FALSE}.
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
\deqn{(k+1)/(n+1) ( (x/X_{n-k,n})^{-1/\gamma_k} -  R_{r,k,n}^{1/\hat{\gamma}_k} ) / (1- R_{r,k,n}^{1/\hat{\gamma}_k})}
with \eqn{R_{r,k,n} = X_{n-k,n} / X_{n-r+1,n}} and \eqn{\hat{\gamma}_k} the Hill estimates adapted for truncation.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates.}
\item{q}{The used large quantile.}
}
\references{
Beirlant, J., Fraga Alves, M.I., Gomes, M.I. and Meerschaert, M.M. (2015). "Extreme Value Statistics for Truncated
Pareto-type Distributions", \emph{http://arxiv.org/abs/1410.4097}.
}
\author{
Tom Reynkens based on \code{R} code of Dries Cornilly.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trHill}}, \code{\link{trQuant}}, \code{\link{Prob}}, \code{\link{trcProb}}, 
\code{\link{trciProb}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000,shape=shape,endpoint=qpareto(0.99,shape=shape))

# Truncated Hill estimator
trh <- trHill(X,plot=TRUE,ylim=c(0,2))

# Small probability
trProb(X,gamma=trh$gamma,q=8,plot=TRUE)
}
