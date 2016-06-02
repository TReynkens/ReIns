\name{Scale}
\alias{Scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale estimator
}
\description{
Computes the estimator for the scale parameter as described in Beirlant et al. (201x).
}
\usage{
Scale(data, gamma = NULL, logk = FALSE, plot = FALSE, add = FALSE, 
      main = "Estimates of scale parameter", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI. When \code{NULL} (the default value), Hill estimates are computed.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of scale parameter"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The scale estimates are computed based on the following model for the CDF:
\eqn{1-F(x) = A x^{-1/\gamma}}, where \eqn{A:= C^{1/\gamma}} is the scale parameter:
\deqn{\hat{A}_{k,n}=(k+1)/(n+1) X_{n-k,n}^(1/H_{k,n})}
where \eqn{H_{k,n}} are the Hill estimates.

See Section 4.2.1 of Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{A}{Vector of the corresponding scale estimates.}
\item{C}{Vector of the corresponding estimates for \eqn{C}, see details.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Schoutens, W., De Spiegeleer, J., Reynkens, T. and Herrmann, K. (201x). 
"Hunting for Black Swans in the European Banking Sector." \emph{}, , --.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ScaleEPD}}, \code{\link{Scale.2o}}, \code{\link{Hill}}
}
\examples{
data(secura)

# Hill estimator
H <- Hill(secura$size)

# Scale estimator
S <- Scale(secura$size, gamma=H$gamma, plot=FALSE)

# Plot logarithm of scale          
plot(S$k,log(S$A), xlab="k", ylab="log(Scale)", type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

