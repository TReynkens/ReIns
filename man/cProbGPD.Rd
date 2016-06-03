\name{cProbGPD}
\alias{cProbGPD}
\alias{cReturnGPD}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using censored GPD-MLE
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the GPD-ML estimator adapted for right censoring.
}
\usage{
cProbGPD(data, censored, gamma1, sigma1, q, plot = FALSE, add = FALSE,
        main = "Estimates of small exceedance probability", ...)

cReturnGPD(data, censored, gamma1, sigma1, q, plot = FALSE, add = FALSE,
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
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{cGPDmle}}.
}  
  \item{sigma1}{
Vector of \eqn{n-1} estimates for \eqn{\sigma_1} obtained from \code{\link{cGPDmle}}.
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{cProbGPD}
and \code{"Estimates of return period"} for \code{cReturnGPD}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The probability is estimated as
\deqn{ \hat{P}(X>q)=(1-km) * (1+ \hat{\gamma}_1/a_{k,n} * (q-Z_{n-k,n}))^{-1/\hat{\gamma}_1}} with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{\hat{\gamma}_1} the generalised Hill estimator adapted for right censoring and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}. The value \eqn{a} is defined as
\deqn{a_{k,n} = \hat{\sigma}_1 / \hat{p}_k} with \eqn{\hat{\sigma}_1} the ML estimate for \eqn{\sigma_1}
and \eqn{\hat{p}_k} the proportion of the \eqn{k} largest observations that is non-censored.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{cProbGPD}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{cReturnGPD}.}
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
\code{\link{cQuantGPD}}, \code{\link{cGPDmle}}, \code{\link{ProbGPD}}, \code{\link{Prob}}, \code{\link{KaplanMeier}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1, length(secura$size), replace=TRUE)

# GPD-MLE estimator adapted for right censoring
cpot <- cGPDmle(secura$size, censored=censored, plot=TRUE)

# Exceedance probability
q <- 10^7
cProbGPD(secura$size, gamma1=cpot$gamma1, sigma1=cpot$sigma1,
         censored=censored, q=q, plot=TRUE)
         
# Return period
cReturnGPD(secura$size, gamma1=cpot$gamma1, sigma1=cpot$sigma1,
         censored=censored, q=q, plot=TRUE)         
}
