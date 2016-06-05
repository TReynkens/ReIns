\name{ProbReg}
\alias{ProbReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small tail probability in regression
}
\description{
Estimator of small tail probability \eqn{1-F_i(q)} in the regression case where \eqn{\gamma} is constant and the regression modelling is thus only solely placed on the scale parameter.
}
\usage{
ProbReg(Z, A, q, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Vector of \eqn{n} observations (from the response variable). 
}
  \item{A}{
Vector of \eqn{n-1} estimates for \eqn{A(i/n)} obtained from \code{\link{ScaleReg}}.
}
  \item{q}{
The used large quantile (we estimate \eqn{P(X_i>q)}) for \eqn{q} large).
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The estimator is defined as
\deqn{1-\hat{F}_i(q) = \hat{A}(i/n) (k+1)/(n+1) (q/Z_{n-k,n})^{-1/H_{k,n}},}
with \eqn{H_{k,n}} the Hill estimator. Here, it is assumed that we have equidistant covariates \eqn{x_i=i/n}.

See Section 4.4 in Albrecher et al. (2016) for more details. 
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{QuantReg}}, \code{\link{ScaleReg}}, \code{\link{Prob}}
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

