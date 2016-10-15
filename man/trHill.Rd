\name{trHill}
\alias{trHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill estimator for upper truncated data
}
\description{
Computes the Hill estimator for positive extreme value indices, adapted for upper truncation, as a function of the tail parameter \eqn{k} (Aban et al. 2006; Beirlant et al., 2016). Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
trHill(data, r = 1, tol = 1e-08, maxiter = 100, logk = FALSE,
       plot = FALSE, add = FALSE, main = "Estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
}
  \item{tol}{
Numerical tolerance for stopping criterion used in Newton-Raphson iterations, default is \code{1e-08}.
}
  \item{maxiter}{
Maximum number of Newton-Raphson iterations, default is \code{100}.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The truncated Hill estimator is the MLE for \eqn{\gamma} under the truncated Pareto distribution.

To estimate the EVI using the truncated Hill estimator an equation needs to be solved. Beirlant et al. (2016) propose
to use Newton-Raphson iterations to solve this equation. We take the trimmed Hill estimates as starting values for this algorithm. The trimmed Hill estimator is defined as
\deqn{H_{r,k,n} = 1/(k-r+1) \sum_{j=r}^k \log(X_{n-j+1,n})-\log(X_{n-k,n})   for 1\le r < k < n}
and is a basic extension of the Hill estimator for upper truncated data (the ordinary Hill estimator is obtained for \eqn{r=1}).

The equation that needs to be solved is
\deqn{H_{r,k,n} = \gamma + R_{r,k,n}^{1/\gamma} \log(R_{r,k,n}) / (1-R_{r,k,n}^{1/\gamma})}
with \eqn{R_{r,k,n} = X_{n-k,n} / X_{n-r+1,n}}.

See Beirlant et al. (2016) or Section 4.2.3 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding estimates for \eqn{\gamma}.}
\item{H}{Vector of corresponding trimmed Hill estimates.}
}
\references{
Aban, I.B., Meerschaert, M.M. and Panorska, A.K. (2006). "Parameter Estimation for the Truncated
Pareto Distribution." \emph{Journal of the American Statistical Association}, 101, 270--277.

Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Fraga Alves, M.I. and Gomes, M.I. (2016). "Tail fitting for Truncated and Non-truncated Pareto-type Distributions." \emph{Extremes}, 19, 429--462.
}
\author{
Tom Reynkens based on \code{R} code of Dries Cornilly.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{trDT}}, \code{\link{trEndpoint}}, \code{\link{trProb}}, \code{\link{trQuant}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000, shape=shape, endpoint=qpareto(0.99, shape=shape))

# Truncated Hill estimator
trh <- trHill(X, plot=TRUE, ylim=c(0,2))
}


