\name{ciLognormalQQ}
\alias{ciLognormalQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-normal QQ-plot for interval censored data
}
\description{
Log-normal QQ-plot adapted for interval censored data.
}
\usage{
ciLognormalQQ(Z, I, censored, plot = TRUE, main = "Log-normal QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals. In an 
insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a log-normal QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Log-normal QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The log-normal QQ-plot adapted for interval censoring is given by
\deqn{( \Phi^{-1}(F_{TB}(Z_{j,n})),  log(Z_{j,n}) )}
for \eqn{j=1,...,n-1,}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{\Phi^{-1}} the quantile function of the standard normal distribution and \eqn{F_{TB}} the Turnbull estimator for the CDF.
Hence, it has the same empirical quantiles as an ordinary log-normal QQ-plot but replaces the theoretical quantiles \eqn{ \Phi^{-1}(j/(n+1))} by \eqn{\Phi^{-1}(F_{TB}(Z_{j,n}))}.
}
\value{
A list with following components:
\item{lqq.the}{Vector of the theoretical quantiles, see Details.}
\item{lqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LognormalQQ}}, \code{\link{cLognormalQQ}}, \code{\link{ciExpQQ}}, \code{\link{ciParetoQQ}}, \code{\link{Turnbull}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Incurred values
I <- secura$size+rchisq(length(secura$size),2)

# Log-normal QQ-plot adapted for interval censoring
ciLognormalQQ(secura$size,I,censored=censored)
}