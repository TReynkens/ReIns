\name{Turnbull}
\alias{Turnbull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turnbull estimator
}
\description{
Computes the Turnbull estimator for the CDF of interval censored data.
}
\usage{
Turnbull(x, L, R, censored)
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
}
\details{
We consider the random right censoring model where one observes \eqn{L \le R}
and where the variable of interest \eqn{X} lies between \eqn{L} and \eqn{R}.

This function calls \code{survfit.formula} from \pkg{survival}.
}
\value{
A vector of length \code{length(x)} containing the Turnbull estimator evaluated in the elements of \code{x}.
}
\references{
Turnbull, B. W. (1974). "Nonparametric Estimation of a Survivorship Function with Doubly Censored Data." \emph{Journal of the American Statistical Association}, 69, 169--173.
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
R <- c(1,2.5,3,4,5.5,6,7.5,8.25,9,10.5)
censored <- c(0,1,0,0,1,0,1,1,0,1)

x <- seq(0,12,0.1)

# Turnbull estimator
plot(x,Turnbull(x, L, R, censored),type="s", ylab="Turnbull estimator")
}
