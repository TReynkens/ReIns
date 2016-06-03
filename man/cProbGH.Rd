\name{cProbGH}
\alias{cProbGH}
\alias{cReturnGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using censored generalised Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the estimates for the EVI obtained from the generalised Hill estimator adapted for right censoring.
}
\usage{
cProbGH(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)

cReturnGH(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
        main = "Estimates of return period", ...)        
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
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{cgenHill}}.
}
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} or \eqn{1/P(X>q)} for \eqn{q} large).
}
 \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{cProbGH}
and \code{"Estimates of return period"} for \code{cReturnGH}.

}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The probability is estimated as
\deqn{ \hat{P}(X>q)=(1-km) * (1+ \hat{\gamma}_1/a_{k,n} * (q-Z_{n-k,n}))^{-1/\hat{\gamma}_1}} with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{\hat{\gamma}_1} the generalised Hill estimator adapted for right censoring and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}. The value \eqn{a} is defined as
\deqn{a_{k,n} = Z_{n-k,n} H_{k,n} (1-\min(\hat{\gamma}_1,0)) / \hat{p}_k} with \eqn{H_{k,n}} the ordinary Hill estimator
and \eqn{\hat{p}_k} the proportion of the \eqn{k} largest observations that is non-censored.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{cProbGH}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{cReturnGH}.}
\item{q}{The used large quantile.}
}
\references{
Einmahl, J.H.J., Fils-Villetard, A. and Guillou, A. (2008). "Statistics of Extremes Under Random Censoring."
\emph{Bernoulli}, 14, 207--227.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cQuantGH}}, \code{\link{cgenHill}}, \code{\link{ProbGH}}, \code{\link{cProbMOM}}, \code{\link{KaplanMeier}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1, length(secura$size), replace=TRUE)

# Generalised Hill estimator adapted for right censoring
cghill <- cgenHill(secura$size, censored=censored, plot=TRUE)

# Small exceedance probability
q <- 10^7
cProbGH(secura$size, censored=censored, gamma1=cghill$gamma1, q=q, plot=TRUE)

# Return period
cReturnGH(secura$size, censored=censored, gamma1=cghill$gamma1, q=q, plot=TRUE)
}