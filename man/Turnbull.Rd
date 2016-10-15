\name{Turnbull}
\alias{Turnbull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turnbull estimator
}
\description{
Computes the Turnbull estimator for the survival function of interval censored data.
}
\usage{
Turnbull(x, L, R, censored, trunclower = 0, truncupper = Inf,
         conf.type = "plain", conf.int = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector with points to evaluate the estimator in.
}
  \item{L}{
Vector of length \eqn{n} with the left boundary of the intervals.
}
  \item{R}{
Vector of length \eqn{n} with the right boundary of the intervals.
}
  \item{censored}{
Vector of \eqn{n} logicals indicating if an observation is interval censored. 
}
\item{trunclower}{
Lower truncation point, default is 0.
}
\item{truncupper}{
Upper truncation point, default is \code{Inf}.
}
  \item{conf.type}{
Type of confidence interval, see \code{\link[survival]{survfit.formula}}. Default is \code{"plain"}.
}
  \item{conf.int}{
Confidence level of the two-sided confidence interval, see \code{\link[survival]{survfit.formula}}. Default is \code{0.95}.
}
}
\details{
We consider the random interval censoring model where one observes \eqn{L \le R}
and where the variable of interest \eqn{X} lies between \eqn{L} and \eqn{R}.

Right censored data should be entered as \code{L=l} and \code{R=truncupper}, and right censored data should be entered as \code{L=trunclower} and \code{R=r}.

This function calls \code{survfit.formula} from \pkg{survival}.

See Section 4.3.2 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{surv}{A vector of length \code{length(x)} containing the Turnbull estimator evaluated in the elements of \code{x}.}
\item{fit}{The output from the call to \code{\link[survival]{survfit.formula}}, an object of class \code{survfit}.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Turnbull, B. W. (1974). "Nonparametric Estimation of a Survivorship Function with Doubly Censored Data." \emph{Journal of the American Statistical Association}, 69, 169--173.

Turnbull, B. W. (1976). "The Empirical Distribution Function with Arbitrarily
Grouped, Censored and Truncated Data." \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, 38, 290--295.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[survival]{survfit.formula}}, \code{\link{KaplanMeier}}
}
\examples{
L <- 1:10
R <- c(1, 2.5, 3, 4, 5.5, 6, 7.5, 8.25, 9, 10.5)
censored <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1)

x <- seq(0, 12, 0.1)

# Turnbull estimator
plot(x, Turnbull(x, L, R, censored)$cdf, type="s", ylab="Turnbull estimator")
}
