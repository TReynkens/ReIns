\name{cMoment}
\alias{cMoment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MOM estimator for right censored data
}
\description{
Computes the Method of Moment estimates adapted for right censored data.
}
\usage{
cMoment(data, censored, logk = FALSE, plot = FALSE, add = FALSE, 
        main = "Moment estimates of the EVI", ...)
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
Logical indicating if the estimates are plotted as a function of \eqn{\log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma_1} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma_1} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Moment estimates of the EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The moment estimator adapted for right censored data is equal to the ordinary moment estimator divided by the proportion of the \eqn{k} largest observations that is non-censored.

This estimator is only suitable for right censored data.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma1}{Vector of the corresponding moment estimates.}
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
\code{\link{Moment}}, \code{\link{cProbMOM}}, \code{\link{cQuantMOM}}
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

# Moment estimator adapted for right censoring
cmom <- cMoment(Z, censored=censored, plot=TRUE)
}
