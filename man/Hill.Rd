\name{Hill}
\alias{Hill}

\title{
Hill estimator
}
\description{
Computes the Hill estimator for positive extreme value indices (Hill, 1975) as a function of the tail parameter \eqn{k}.
Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
Hill(data, k = TRUE, plot = FALSE, add = FALSE, main = "Hill estimates of EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations.  
}
  \item{k}{
Logical indicating if the Hill estimates are plotted as a function of the tail parameter \eqn{k} (\code{k=TRUE}) or as a function of \eqn{log X_{n-k}}. Default is \code{TRUE}.
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Hill estimates of EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}

\details{
The Hill estimator can be seen as the estimator of slope in the upper right corner (\eqn{k} last points) of the Pareto QQ-plot when using constrained least squares (the regression line has to pass through the point \eqn{-\log((k+1)/(n+1)),\log X_{n-k}}). See Section 4.2 of Beirlant et al. (2004) for more details.
}

\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding Hill estimates.}
}

\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley.

Hill, B. M. (1975). "A simple general approach to inference about the tail of a
distribution." \emph{Annals of Statistics}, 3, 1163--1173.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ParetoQQ}}, \code{\link{Hill.2oQV}}, \code{\link{genHill}}
}
\examples{
data(norwegianfire)

# Plot Hill estimates as a function of k
Hill(norwegianfire$size[norwegianfire$year==76],plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
