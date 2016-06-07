\name{icHill}
\alias{icHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill estimator for interval censored data
}
\description{
Computes the Hill estimator for positive extreme value indices, adapted for interval censoring, as a function of the tail parameter \eqn{k}. Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
icHill(L, U, censored, trunclower = 0, truncupper = Inf, 
       logk = FALSE, plot = TRUE, add = FALSE, main = "Hill estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
Vector of length \eqn{n} with the lower boundary of the intervals for interval censored data or the observed data for right censored data.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundary of the intervals.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{trunclower}{
Lower truncation point. Default is 0.
}
  \item{truncupper}{
Upper truncation point. Default is \code{Inf} (no upper truncation).
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
Title for the plot, default is \code{"Hill estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
This estimator is given by 
\deqn{H^{TB}(x)=(\int_x^{\infty} (1-\hat{F}^{TB}(u))/u du)/(1-\hat{F}^{TB}(x)),} 
where \eqn{\hat{F}^{TB}} is the Turnbull estimator for the CDF.
More specifically, we use the values \eqn{x=\hat{Q}^{TB}(p)} for \eqn{p=1/(n+1), \ldots, (n-1)/(n+1)} where
\eqn{\hat{Q}^{TB}(p)} is the empirical quantile function corresponding to the Turnbull estimator.
We then denote \deqn{H^{TB}_{k,n}=H^{TB}(x_{n-k,n})} with \deqn{x_{n-k,n}=\hat{Q}^{TB}((n-k)/(n+1))=\hat{Q}^{TB}(1-(k+1)/(n+1)).}

Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and left censored data should be entered as \code{L=trunclower} and \code{U=u}. 

If the \pkg{interval} package is installed, the \code{\link[interval]{icfit}} function is used to compute the Turnbull estimator. Otherwise, \code{\link[survival]{survfit.formula}} from \pkg{survival} is used.

Use \code{\link{Hill}} for non-censored data or \code{\link{cHill}} for right censored data. 

See Section 4.3 in Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding Hill estimates.}
\item{X}{Vector of thresholds \eqn{x_{n-k,n}} used when estimating \eqn{\gamma}.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cHill}}, \code{\link{Hill}}, \code{\link{MeanExcess_TB}}, \code{\link{icParetoQQ}}, \code{\link{Turnbull}}, \code{\link[interval]{icfit}}
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

# Hill estimator adapted for interval censoring
icHill(Z, U, censored, ylim=c(0,1))

# Hill estimator adapted for right censoring
cHill(Z, censored, lty=2, add=TRUE)

# True value of gamma
abline(h=1/2, lty=3, col="blue")

# Legend
legend("topright", c("icHill", "cHill"), lty=1:2)
}
