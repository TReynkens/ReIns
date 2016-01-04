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
Scale(data, gamma = NULL, plot = FALSE, add = FALSE, 
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
\eqn{1-F(x) = A x^{-1/\gamma}}, where \eqn{A:= C^{1/\gamma}} is the scale parameter.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{A}{Vector of the corresponding scale estimates.}
\item{C}{Vector of the corresponding estimates for \eqn{C}, see details.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Beirlant, J., De Spiegeleer, J., Herrmann, K., Reynkens, T. and Schoutens, W. (201x). 
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

