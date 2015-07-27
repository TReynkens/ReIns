\name{trciHill}
\alias{trciHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill estimator for interval censored upper truncated data
}
\description{
Computes the Hill estimator for positive extreme value indices, adapted for intervl censoring and upper truncation, as a function of the tail parameter. Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
trciHill(Z, I, censored, r = 1, threshold = NULL, lower = 0.01, upper = 2, 
        up = 0.2, maxup = 50, warnings = FALSE, plot = FALSE, add = FALSE, 
        main = "Estimates for the EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals. In an 
insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values.
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
 \item{warnings}{
Logical indicating if possible warnings from \code{\link[stats]{uniroot}} are shown, default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"EPD estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
To estimate the EVI \eqn{\gamma_1}, the root of following equation has to be computed:
\deqn{H_{t,n}^c = \gamma_1 + g(\gamma_1;Z,I,\delta),}
where \eqn{t} is the threshold, \eqn{H_{t,n}^c} the Hill estimator adapted for right censoring, \eqn{\delta} the non-censoring indicator (the complement of \code{censoring}) and \eqn{g} 
a function that depends on \eqn{gamma_1}, see Section 3.2.5 of Beirlant et al. (201x) for more details.

We compute this root using \code{\link[stats]{uniroot}}. Therefore, the arguments \code{lower} and \code{upper} need to be specified. If the function (of which we compute the roots) is equal in \code{lower} and \code{upper}, we increase \code{upper}
by \code{up} until the function values in \code{lower} and \code{upper} are different or until we increased it \code{maxup} times. In the latter case, we return \code{NA}.


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
\code{\link{trciProb}}, \code{\link{ciHill}}, \code{\link{trHill}}, \code{\link{trcHill}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=100,shape=shape,endpoint=qpareto(0.99,shape=shape))

# Incurred values
I <- X+rchisq(length(X),2)

# Truncated Hill estimator
trh <- trHill(X,plot=TRUE,ylim=c(0,2))

# Censoring indicator
censored <- sample(0:1,length(X),replace=TRUE)

# Hill estimator adapted for right censoring and 
# upper truncation
trcih <- trciHill(X,I,censored=censored,plot=TRUE,ylim=c(0,2))
}