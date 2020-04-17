\name{cProb}
\alias{cProb}
\alias{cReturn}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and large return periods using censored Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or large return period \eqn{1/P(X>q)} using the estimates for the EVI obtained from the Hill estimator adapted for right censoring.
}
\usage{
cProb(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
      main = "Estimates of small exceedance probability", ...)
      
cReturn(data, censored, gamma1, q, plot = FALSE, add = FALSE, 
        main = "Estimates of large return period", ...)      
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{cProb}
and \code{"Estimates of large return period"} for \code{cReturn}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The probability is estimated as
\deqn{ \hat{P}(X>q)=(1-km) \times (q/Z_{n-k,n})^{-1/H_{k,n}^c}} with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{H_{k,n}^c}
the Hill estimator adapted for right censoring and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{cProb}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{cReturn}.}
\item{q}{The used large quantile.}
}
\references{
 Beirlant, J., Guillou, A., Dierckx, G. and Fils-Villetard, A. (2007). "Estimation of the Extreme Value Index and Extreme Quantiles Under Random Censoring." \emph{Extremes}, 10, 151--174.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cHill}}, \code{\link{cQuant}}, \code{\link{Prob}}, \code{\link{KaplanMeier}}
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

# Small exceedance probability
q <- 10
cProb(Z, censored=censored, gamma1=chill$gamma1, q=q, plot=TRUE)

# Return period
cReturn(Z, censored=censored, gamma1=chill$gamma1, q=q, plot=TRUE)
}

