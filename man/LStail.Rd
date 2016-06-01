
\name{LStail}
\alias{LStail}
\alias{TSfraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least Squares tail estimator
}
\description{
Computes the Least Squares (LS) estimates of the EVI based on the last \eqn{k} observations of the generalised QQ-plot.
}
\usage{
LStail(data, rho = -1, lambda = 0.5, logk = FALSE, plot = FALSE, add = FALSE, 
            main = "LS estimates of EVI", ...)
            
TSfraction(data, rho = -1, lambda = 0.5, logk = FALSE, plot = FALSE, add = FALSE, 
            main = "LS estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{rho}{
Estimate for \eqn{\rho}, or \code{NULL} when \eqn{\rho} needs to be estimated using the method of Beirlant et al. (2002). Default is \code{-1}. 
}
  \item{lambda}{
Parameter used in the method of Beirlant et al. (2002), only used when \code{rho=NULL}. Default is \code{0.5}.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of the \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
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
We estimate \eqn{\gamma} (EVI) and \eqn{b} using least squares on the following regression model (Beirlant et al., 2005): \eqn{Z_j = \gamma + b(n/k) (j/k)^{-\rho} + \epsilon_j} with \eqn{Z_j = (j+1) \log(UH_{j,n}/UH_{j+1,n})} and \eqn{UH_{j,n}=X_{n-j,n}H_{j,n}} where \eqn{H_{j,n}} is the Hill estimator with threshold \eqn{X_{n-j,n}}.
See Section 5.8 of Beirlant et al. (2004) for more details.

The function \code{TSfraction} is included for compatibility with the old \code{S-Plus} code.

}
\value{
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding LS estimates for the EVI.}
\item{b}{Vector of the corresponding LS estimates for b.}
\item{rho}{Vector of the estimates for \eqn{\rho} when \code{rho=NULL} or the given input for \code{rho} otherwise.}
}
\references{

Beirlant, J., Dierckx, G. and Guillou, A. (2005). "Estimation of the Extreme Value Index and Regression on Generalized Quantile Plots." \emph{Bernoulli}, 11, 949--970.

Beirlant, J., Dierckx, G., Guillou, A. and Starica, C. (2002). "On Exponential Representations of Log-spacing of Extreme Order Statistics." \emph{Extremes}, 5, 157--180.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.


}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{genQQ}}
}
\examples{
data(soa)

# LS tail estimator
LStail(soa$size,plot=TRUE,ylim=c(0,0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

