\name{GPDmle}
\alias{GPDmle}
\alias{POT}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GPD-ML estimator
}
\description{
Fit the Generalised Pareto Distribution (GPD) to the exceedances (peaks) over a threshold using Maximum Likelihood Estimation (MLE). Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
GPDmle(data, start = c(0.1,1), warnings = FALSE, logk = FALSE, 
       plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...)

POT(data, start = c(0.1,1), warnings = FALSE, logk = FALSE, 
    plot = FALSE, add = FALSE, main = "POT estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{start}{
Vector of length 2 containing the starting values for the optimisation. The first element
is the starting value for the estimator of \eqn{\gamma} and the second element is the starting value for the estimator of \eqn{\sigma}. Default is \code{c(0.1,1)}.
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
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
Title for the plot, default is \code{"POT estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The \code{POT} function is the same function but with a different name for compatibility with the old \code{S-Plus} code. 

For each value of \code{k}, we look at the exceedances over the \eqn{k+1}th largest observations:
\eqn{X_{n-k+j,n}-X_{n-k,n}} for \eqn{j=1,...,k}, with \eqn{X_{j,n}} the \eqn{j}th largest observation and \eqn{n} the sample size. The GPD is then fitted to these k exceedances using MLE which yields estimates for the parameters of the GPD: \eqn{\gamma} and \eqn{\sigma}. See Section 5.3 of Beirlant et al. (2004) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding MLE estimates for the \eqn{\gamma} parameter of the GPD.}
\item{sigma}{Vector of the corresponding MLE estimates for the \eqn{\sigma} parameter of the GPD.}
}
\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur and \code{R} code from Klaus Herrmann.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GPDfit}}, \code{\link{GPDresiduals}}, \code{\link{EPD}}
}

\examples{
data(soa)

# Look at last 500 obs of SOA data
SOAdata <- sort(soa$size)[length(soa$size)-(0:499)]

# Plot GPD-ML estimates as a function of k
GPDmle(SOAdata, plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

