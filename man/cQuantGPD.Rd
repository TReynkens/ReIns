\name{cQuantGPD}
\alias{cQuantGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of large quantiles using censored GPD-MLE
}
\description{
Computes estimates of large quantiles \eqn{Q(1-p)} using the estimates for the EVI obtained from the GPD-ML estimator adapted for right censoring.
}
\usage{
cQuantGPD(data, censored, gamma1, sigma1, p, plot = FALSE, add = FALSE, 
          main = "Estimates of extreme quantile", ...)
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
Title for the plot, default is \code{"Estimates of extreme quantile"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The quantile is estimated as 
\deqn{\hat{Q}(1-p)= Z_{n-k,n} + a_{k,n} ( ( (1-km)/p)^{\hat{\gamma}_1} -1 ) / \hat{\gamma}_1)} 
ith \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{\hat{\gamma}_1} the generalised Hill estimator adapted for right censoring and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}. The value \eqn{a} is defined as
\deqn{a_{k,n} = \hat{\sigma}_1 / \hat{p}_k} with \eqn{\hat{\sigma}_1} the ML estimate for \eqn{\sigma_1}
and \eqn{\hat{p}_k} the proportion of the \eqn{k} largest observations that is non-censored.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
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
\code{\link{cProbGPD}}, \code{\link{cGPDmle}}, \code{\link{QuantGPD}}, \code{\link{Quant}}, \code{\link{KaplanMeier}}
}
\examples{
# Set seed
set.seed(29072016)

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X, Y)

# Censoring indicator
censored <- (X>Y)

# GPD-MLE estimator adapted for right censoring
cpot <- cGPDmle(Z, censored=censored, plot=TRUE)

# Large quantile
p <- 10^(-4)
cQuantGPD(Z, gamma1=cpot$gamma1, sigma1=cpot$sigma1,
         censored=censored, p=p, plot=TRUE)
}         