\name{Hill.2oQV}
\alias{Hill.2oQV}

\title{
Bias-reduced MLE (Quantile view)
}
\description{
Computes bias-reduced ML estimates of gamma based on the quantile view.
}
\usage{
Hill.2oQV(data, start = c(1,1,1), warnings = FALSE, logk = FALSE, 
          plot = FALSE, add = FALSE, main = "Estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{start}{
A vector of length 3 containing starting values for the first numerical optimisation (see Details). The elements
are the starting values for the estimators of \eqn{\gamma}, \eqn{\mu} and \eqn{\sigma}, respectively. Default is \code{c(1,1,1)}.
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
  \item{logk}{
Logical indicating if the estimates are plotted as a function of \eqn{\log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.1 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the ML estimates for the EVI for each value of \eqn{k}.}
\item{b}{Vector of the ML estimates for the parameter \eqn{b} in the regression model for each value of \eqn{k}.}
\item{beta}{Vector of the ML estimates for the parameter \eqn{\beta} in the regression model for each value of \eqn{k}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Dierckx, G., Goegebeur Y. and Matthys, G. (1999). "Tail Index Estimation and an Exponential Regression Model." \emph{Extremes}, 2, 177--200.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur and \code{R} code from Klaus Herrmann.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
data(norwegianfire)

# Plot bias-reduced MLE (QV) as a function of k
Hill.2oQV(norwegianfire$size[norwegianfire$year==76],plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

