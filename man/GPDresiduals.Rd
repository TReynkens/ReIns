\name{GPDresiduals}
\alias{GPDresiduals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GPD residual plot
}
\description{
Residual plot to check GPD fit for peaks over a threshold.
}
\usage{
GPDresiduals(data, t, gamma, sigma, plot = TRUE, 
             main = "GPD residual plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{t}{
The used threshold.
}
  \item{gamma}{
Estimate for the EVI obtained from \code{GPDmle}.
}
  \item{sigma}{
Estimate for \eqn{\sigma} obtained from \code{GPDmle}.
}
  \item{plot}{
Logical indicating if the residuals should be plotted, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"GPD residual plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Consider the POT values \eqn{Y=X-t} and the transformed variable
\deqn{R= 1/\gamma \log(1+\gamma/\sigma Y), } when \eqn{\gamma \neq 0} and
\deqn{R = Y/\sigma,}
otherwise. 
We can assess the goodness-of-fit of the GPD when modelling POT values \eqn{Y=X-t} by 
constructing an exponential QQplot of the transformed variable \eqn{R} since \eqn{R} is standard exponentially distributed if \eqn{Y} follows the GPD.
}
\value{
A list with following components:
\item{res.the}{Vector of the theoretical quantiles from a standard exponential distribution.}
\item{res.emp}{Vector of the empirical quantiles of \eqn{R}, see Details.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GPDfit}}, \code{\link{ExpQQ}}
}
\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Plot POT-MLE estimates as a function of k
pot <- GPDmle(SOAdata, plot=TRUE)

# Residual plot
}
