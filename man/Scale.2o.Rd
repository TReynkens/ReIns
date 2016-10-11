\name{Scale.2o}
\alias{Scale.2o}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bias-reduced scale estimator using second order Hill estimator
}
\description{
Computes the bias-reduced estimator for the scale parameter using the second order Hill estimator.
}
\usage{
Scale.2o(data, gamma, b, beta, logk = FALSE, plot = FALSE, add = FALSE, 
        main = "Estimates of scale parameter", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{Hill.2oQV}.
}
\item{b}{
Vector of \eqn{n-1} estimates for \eqn{B} obtained from \code{Hill.2oQV}.
}
  \item{beta}{
Vector of \eqn{n-1} estimates for \eqn{\beta} obtained from \code{Hill.2oQV}.
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
\eqn{1-F(x) = A x^{-1/\gamma} ( 1+ bx^{-\beta}(1+o(1)) )}, where \eqn{A:= C^{1/\gamma}} is the scale parameter.

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

Beirlant, J., De Spiegeleer, J., Schoutens, W., Reynkens, T. and Herrmann, K. (2016).
"Hunting for Black Swans in the European Banking Sector using Extreme Value Analysis."
In: Jan Kallsen and Antonis Papapantoleon (eds.), \emph{Advanced Modeling in Mathematical
Finance: In honour of Ernst Eberlein}, Springer International Publishing, Switzerland,
accepted.
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{Scale}}, \code{\link{ScaleEPD}}, \code{\link{Hill.2oQV}}
}
\examples{
data(secura)

# Hill estimator
H <- Hill(secura$size)
# Bias-reduced Hill estimator
H2o <- Hill.2oQV(secura$size)

# Scale estimator
S <- Scale(secura$size, gamma=H$gamma, plot=FALSE)
# Bias-reduced scale estimator
S2o <- Scale.2o(secura$size, gamma=H2o$gamma, b=H2o$b, 
          beta=H2o$beta, plot=FALSE)

# Plot logarithm of scale             
plot(S$k,log(S$A), xlab="k", ylab="log(Scale)", type="l")
lines(S2o$k,log(S2o$A), lty=2)
}
