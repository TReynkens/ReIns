\name{ScaleReg}
\alias{ScaleReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scale estimator in regression
}
\description{
Estimator of the scale parameter in the regression case where \eqn{\gamma} is constant and the regression modelling is thus placed solely on the scale parameter.
}
\usage{
ScaleReg(s, Z, kernel = c("normal", "uniform", "triangular", "epanechnikov", "biweight"), 
         h, plot = TRUE, add = FALSE, main = "Estimates of scale parameter", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
Point to evaluate the scale estimator in.
}
  \item{Z}{
Vector of \eqn{n} observations (from the response variable). 
}
  \item{kernel}{
The kernel used in the estimator. One of \code{"normal", "uniform", "triangular", "epanechnikov", "biweight"}.
}
  \item{h}{
The bandwidth used in the kernel function.
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
The scale estimator is computed as
\deqn{\hat{A}(s) = 1/(k+1) \sum_{i=1}^n 1_{Z_i>Z_{n-k,n}} K_h(s-i/n)}
with \eqn{K_h(x)=K(x/h)/h,} \eqn{K} the kernel function and \eqn{h} the bandwidth.
Here, it is assumed that we have equidistant covariates \eqn{x_i=i/n}.

See Section 4.4.1 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{A}{Vector of the corresponding scale estimates.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ProbReg}}, \code{\link{QuantReg}}, \code{\link{scale}}, \code{\link{Hill}}
}
\examples{
data(norwegianfire)

Z <- norwegianfire$size[norwegianfire$year==76]

i <- 100
n <- length(Z)

# Scale estimator in i/n
A <- ScaleReg(i/n, Z, h=0.5, kernel = "epanechnikov")$A

# Small exceedance probability
q <- 10^6
ProbReg(Z, A, q, plot=TRUE)

# Large quantile
p <- 10^(-5)
QuantReg(Z, A, p, plot=TRUE)
}
