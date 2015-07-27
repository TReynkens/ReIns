\name{trcHill}
\alias{trcHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill estimator for right censored upper truncated data
}
\description{
Computes the Hill estimator for positive extreme value indices, adapted for right censoring and upper truncation, as a function of the tail parameter. Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
trcHill(data, censored, r = 1, threshold = NULL, lower = 0.01, upper = 2, 
        up = 0.2, maxup = 50, plot = FALSE, add = FALSE,
        main = "Estimates for the EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
}
  \item{threshold}{
A vector of thresholds. If \code{NULL} (default), the order statistics \eqn{Z_{n-k,n}}
for \eqn{k=1,...,n-1} are used.
}
  \item{lower}{
Lower endpoint of the interval to search over, see Details. Default value is \code{0.01}.
}
  \item{upper}{
Upper endpoint of the interval to search over, see Details. Default value is \code{2}.
}
  \item{up}{
Value to increase upper with if \eqn{f(lower)=f(upper)}, see Details. Default is \code{0.02}.
}
  \item{maxup}{
Maximal amount of times \code{upper} is increased, see Details. Default is \code{50}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma_1} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma_1} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
To estimate the EVI \eqn{\gamma_1}, the root of following equation has to be computed:
\deqn{H_{t,n}^c = \gamma_1 + g(\gamma_1;Z,t,\delta),}
where \eqn{t} is the threshold, \eqn{H_{t,n}^c} the Hill estimator adapted for right censoring, \eqn{\delta} the non-censoring indicator (the complement of \code{censoring}) and \eqn{g} 
a function that depends on \eqn{gamma_1}, see Section 3.2.5 of Beirlant et al. (201x) for more details.

We compute this root using \code{\link[stats]{uniroot}}. Therefore, the arguments \code{lower} and \code{upper} need to be specified. If the function (of which we compute the roots) is equal in \code{lower} and \code{upper}, we increase \code{upper}
by \code{up} until the function values in \code{lower} and \code{upper} are different or until we increased it \code{maxup} times. In the latter case, we return \code{NA}.

This function is only suitable for right censored data. For interval censored data, one can use \code{\link{trciHill}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}, only returned when \code{threshold=NULL}.}
\item{threshold}{Vector of the threshold values.}
\item{gamma1}{Vector of the corresponding Hill estimates.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trcProb}}, \code{\link{cHill}}, \code{\link{trHill}}, \code{\link{trciHill}}
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
}
