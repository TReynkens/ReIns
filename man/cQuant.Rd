\name{cQuant}
\alias{cQuant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of large quantiles using censored Hill
}
\description{
Computes estimates of large quantiles \eqn{Q(1-p)} using the estimates for the EVI obtained from the Hill estimator adapted for right censoring.
}
\usage{
cQuant(data, censored, gamma1, p, plot = FALSE, add = FALSE, 
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
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{cHill}}.
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
The quantile is estimated as \deqn{\hat{Q}(1-p)=Z_{n-k,n} \times ( (1-km)/p)^{H_{k,n}^c}} with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{H_{k,n}^c}
the Hill estimator adapted for right censoring and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{Q}{Vector of the corresponding quantile estimates.}
\item{p}{The used exceedance probability.}
}
\references{
 Beirlant, J., Guillou, A., Dierckx, G. and Fils-Villetard, A. (2007). "Estimation of the Extreme Value Index and Extreme Quantiles Under Random Censoring." \emph{Extremes}, 10, 151--174.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cHill}}, \code{\link{cProb}}, \code{\link{Quant}}, \code{\link{KaplanMeier}}
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

# Hill estimator adapted for right censoring
chill <- cHill(Z, censored=censored, plot=TRUE)

# Large quantile
p <- 10^(-4)
cQuant(Z, gamma1=chill$gamma, censored=censored, p=p, plot=TRUE)
}
