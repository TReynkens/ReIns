\name{KaplanMeier}
\alias{KaplanMeier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kaplan-Meier estimator
}
\description{
Computes the Kaplan-Meier estimator for the survival function of right censored data.
}
\usage{
KaplanMeier(x, data, censored, conf.type="plain", conf.int = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector with points to evaluate the estimator in.
}
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{censored}{
Vector of \eqn{n} logicals indicating if an observation is right censored. 
}
  \item{conf.type}{
Type of confidence interval, see \code{\link[survival]{survfit.formula}}. Default is \code{"plain"}.
}
  \item{conf.int}{
Confidence level of the two-sided confidence interval, see \code{\link[survival]{survfit.formula}}. Default is \code{0.95}.
}
}
\details{
We consider the random right censoring model where one observes \eqn{Z = \min(X,C)}
where \eqn{X} is the variable of interest and \eqn{C} is the censoring variable.

This function is merely a wrapper for \code{\link[survival]{survfit.formula}} from \pkg{survival}.

This estimator is only suitable for right censored data. When the data are interval censored, one can use the Turnbull estimator implemented in \code{\link{Turnbull}}.
}
\value{
A list with following components:
\item{surv}{A vector of length \code{length(x)} containing the Kaplan-Meier estimator evaluated in the elements of \code{x}.}
\item{fit}{The output from the call to \code{\link[survival]{survfit.formula}}, an object of class \code{survfit}.}
}
\references{
Kaplan, E. L. and Meier, P. (1958). "Nonparametric Estimation from Incomplete Observations." \emph{Journal of the American Statistical Association}, 53, 457--481.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survival]{survfit.formula}}, \code{\link{Turnbull}}
}
\examples{
data <- c(1, 2.5, 3, 4, 5.5, 6, 7.5, 8.25, 9, 10.5)
censored <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1)

x <- seq(0, 12, 0.1)

# Kaplan-Meier estimator
plot(x, KaplanMeier(x, data, censored)$surv, type="s", ylab="Kaplan-Meier estimator")
}
