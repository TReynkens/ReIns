\name{icParetoQQ}
\alias{icParetoQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pareto quantile plot for interval censored data
}
\description{
Pareto QQ-plot adapted for interval censored data using the Turnbull estimator.
}
\usage{
icParetoQQ(L, U = L, censored, trunclower = 0, truncupper = Inf, 
           plot = TRUE, main = "Pareto QQ-plot", ...)
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
  \item{plot}{
Logical indicating if the quantiles should be plotted in a Pareto QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Pareto QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The Pareto QQ-plot adapted for interval censoring is given by
\deqn{( -\log(1-F^{TB}(x_{j,n})), \log x_{j,n} )}
for \eqn{j=1,\ldots,n-1,} where \eqn{\hat{F}^{TB}} is the Turnbull estimator for the CDF and \eqn{x_{i,n}=\hat{Q}^{TB}(i/(n+1))} with \eqn{\hat{Q}^{TB}(p)} the empirical quantile function corresponding to the Turnbull estimator.

Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and left censored data should be entered as \code{L=trunclower} and \code{U=u}. 

If the \pkg{interval} package is installed, the \code{\link[interval]{icfit}} function is used to compute the Turnbull estimator. Otherwise, \code{\link[survival]{survfit.formula}} from \pkg{survival} is used.

Use \code{\link{ParetoQQ}} for non-censored data or \code{\link{cParetoQQ}} for right censored data. 

See Section 4.3 in Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{pqq.the}{Vector of the theoretical quantiles, see Details.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cParetoQQ}}, \code{\link{ParetoQQ}}, \code{\link{icHill}}, \code{\link{Turnbull}}, \code{\link[interval]{icfit}}
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

# Pareto QQ-plot adapted for interval censoring
icParetoQQ(Z, U, censored)

# Pareto QQ-plot adapted for right censoring
cParetoQQ(Z, censored)
}

