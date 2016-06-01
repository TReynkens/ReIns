\name{MeanExcess_TB}
\alias{MeanExcess_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean excess function using Turnbull estimator
}
\description{
Computes mean excess values using the Turnbull estimator. These mean excess values can then be plotted as a function of the empirical quantiles (computed using the Turnbull estimator) or as a function of the tail parameter \eqn{k}.
}
\usage{
MeanExcess_TB(L, U = L, censored, trunclower = 0, truncupper = Inf, 
              plot = TRUE, k = TRUE, main = "Mean excess plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
Vector of length \eqn{n} with the lower boundary of the intervals for interval censored data or the observed data for right censored data. In an insurance context, these are the paid amounts.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundary of the intervals. In an 
insurance context, these are the incurred values. By default, they are equal to \code{L}.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{trunclower}{
Lower truncation point, default is 0.
}
  \item{truncupper}{
Upper truncation point, default is \code{Inf}.
}
  \item{plot}{
Logical indicating if the mean excess values should be plotted in a mean excess plot, default is \code{TRUE}.
}
  \item{k}{
Logical indicating if the mean excess values are plotted as a function of the tail parameter \eqn{k} (\code{k=TRUE}) or as a function of the empirical quantiles computed using the Turnbull estimator (\code{k=FALSE}). Default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Mean excess plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The mean excess values are given by
\deqn{\hat{e}^{TB}(v)=(\int_v^{\infty} 1-\hat{F}^{TB}(u) du)/(1-\hat{F}^{TB}(u))}
where \eqn{\hat{F}^{TB}} is the Turnbull estimator for the CDF.
More specifically, we use the values \eqn{v=\hat{Q}^{TB}(p)} for \eqn{p=1/(n+1), \ldots, (n-1)/(n+1)} where
\eqn{\hat{Q}^{TB}(p)} is the empirical quantile function corresponding to the Turnbull estimator.
See Albrecher et al. (2016) for more details.

Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and right censored data should be entered as \code{L=trunclower} and \code{U=u}. 

If the \pkg{interval} package is installed, the \code{\link[interval]{icfit}} function is used to compute the Turnbull estimator. Otherwise, \code{\link[survival]{survfit.formula}} from \pkg{survival} is used.

Use \code{\link{MeanExcess}} for non-censored data. 
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \code{k}.}
\item{X}{Vector of the empirical quantiles, computed using the Turnbull estimator, corresponding to \code{(n-k+1)/(n+1)}.}
\item{e}{Vector of the mean excess values corresponding to the tail parameters in \code{k}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{MeanExcess}}
}

\examples{
# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y  <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X,Y)

# Censoring indicator
censored <- (X>Y)

# Right boundary
U <- Z
U[censored] <- Inf

# Mean excess plot
MeanExcess_TB(Z, U, censored, k=FALSE)
}

