\name{cProbEPD}
\alias{cProbEPD}
\alias{cReturnEPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and large return periods using censored EPD
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or large return period \eqn{1/P(X>q)} using the parameters from the EPD fit adapted for right censoring.
}
\usage{
cProbEPD(data, censored, gamma1, kappa1, beta, q, plot = FALSE, add = FALSE,
         main = "Estimates of small exceedance probability", ...)

cReturnEPD(data, censored, gamma1, kappa1, beta, q, plot = FALSE, add = FALSE,
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
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{cEPD}}.
}
  \item{kappa1}{
Vector of \eqn{n-1} estimates for \eqn{\kappa_1} obtained from \code{\link{cEPD}}.
}
  \item{beta}{
Vector of \eqn{n-1} estimates for \eqn{\beta} obtained from \code{\link{cEPD}}.
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
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{cProbEPD}
and \code{"Estimates of large return period"} for \code{cReturnEPD}.

}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The probability is estimated as \deqn{\hat{P}(X>q)= (1-km) \times (1-F(q))} with \eqn{F}
the CDF of the EPD with estimated parameters \eqn{\hat{\gamma}_1}, \eqn{\hat{\kappa}_1} and \eqn{\hat{\tau}=-\hat{\beta}}
and \eqn{km} the Kaplan-Meier estimator for the CDF evaluated in \eqn{Z_{n-k,n}} (the \eqn{(k+1)}-th largest data point).
}

\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{cProbEPD}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{cReturnEPD}.}
\item{q}{The used large quantile.}
}
\references{
Beirlant, J., Bardoutsos, A., de Wet, T. and Gijbels, I. (2016). "Bias Reduced Tail Estimation for Censored Pareto Type Distributions." \emph{Statistics & Probability Letters}, 109, 78--88.
}
\author{
Tom Reynkens.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cEPD}}, \code{\link{ProbEPD}}, \code{\link{Prob}}, \code{\link{KaplanMeier}}
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

# EPD estimator adapted for right censoring
cepd <- cEPD(Z, censored=censored, plot=TRUE)

# Small exceedance probability
q <- 10
cProbEPD(Z, censored=censored, gamma1=cepd$gamma1,
        kappa1=cepd$kappa1, beta=cepd$beta, q=q, plot=TRUE)

# Return period
cReturnEPD(Z, censored=censored, gamma1=cepd$gamma1,
        kappa1=cepd$kappa1, beta=cepd$beta, q=q, plot=TRUE)        
}

