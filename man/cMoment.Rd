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
cMoment(data, censored, plot = FALSE, add = FALSE, main = "Moment estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Moment estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The moment estimator adapted for right censored data is equal to the ordinary moment estimator divided by the proportion of the \eqn{k} largest observations that is non-censored.

This estimator is only suitable for right censored data!
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
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Moment estimator adapted for right censoring
cmom <- cMoment(secura$size,censored=censored,plot=TRUE)
}
