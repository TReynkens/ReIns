\name{cHill}
\alias{cHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill estimator for right censored data
}
\description{
Computes the Hill estimator for positive extreme value indices, adapted for right censoring, as a function of the tail parameter \eqn{k} (Beirlant et al., 2007).
Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
cHill(data, censored, logk = FALSE, plot = FALSE, add = FALSE, 
      main = "Hill estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma_1} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma_1} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Hill estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The Hill estimator adapted for right censored data is equal to the ordinary Hill estimator \eqn{H_{k,n}} divided by the proportion of the \eqn{k} largest observations that is non-censored.

This estimator is only suitable for right censored data, use \code{\link{icHill}} for interval censored data.

See Section 4.4.3 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma1}{Vector of the corresponding Hill estimates.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Guillou, A., Dierckx, G. and Fils-Villetard, A. (2007). "Estimation of the Extreme Value Index and Extreme Quantiles Under Random Censoring." \emph{Extremes}, 10, 151--174.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{icHill}}, \code{\link{cParetoQQ}}, \code{\link{cProb}}, \code{\link{cQuant}}
}
\examples{
# Set seed
set.seed(29072016)

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y  <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X, Y)

# Censoring indicator
censored <- (X>Y)

# Hill estimator adapted for right censoring
chill <- cHill(Z, censored=censored, plot=TRUE)
}
