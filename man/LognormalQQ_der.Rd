\name{LognormalQQ_der}
\alias{LognormalQQ_der}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derivative plot of the log-normal QQ-plot
}
\description{
Computes the derivative plot of the log-normal QQ-plot. These values can be plotted as a function of the data or as a function of the tail parameter \eqn{k}.
}
\usage{
LognormalQQ_der(data, k = FALSE, plot = TRUE, 
                main = "Derivative plot of log-normal QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the derivative values should be plotted, default is \code{TRUE}.
}
  \item{k}{
Logical indicating if the derivative values are plotted as a function of the tail parameter \eqn{k} (\code{k=TRUE}) or as a function of the logarithm of the data (\code{k=FALSE}). Default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Derivative plot of log-normal QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The derivative plot of a log-normal QQ-plot is \deqn{(k, H_{k,n}/N_{k,n})} or \deqn{(\log X_{n-k,n}, H_{k,n}/N_{k,n})} with \eqn{H_{k,n}} the Hill estimates and \deqn{N_{k,n} = (n+1)/(k+1) \phi(\Phi^{-1}(a)) - \Phi^{-1}(a).} Here is \eqn{a=1-(k+1)/(n+1)}, \eqn{\phi} the normal PDF and \eqn{\Phi} the normal CDF.

See Section 4.1 of Albrecher et al. (2017) for more details.
}

\value{
A list with following components:
\item{xval}{Vector of the x-values of the plot (\eqn{k} or \eqn{\log X_{n-k,n}}).}
\item{yval}{Vector of the derivative values.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LognormalQQ}}, \code{\link{Hill}}, \code{\link{MeanExcess}}, \code{\link{ParetoQQ_der}}, \code{\link{WeibullQQ_der}} 
}
\examples{
data(norwegianfire)

# Log-normal QQ-plot for Norwegian Fire Insurance data for claims in 1976.
LognormalQQ(norwegianfire$size[norwegianfire$year==76])

# Derivate plot
LognormalQQ_der(norwegianfire$size[norwegianfire$year==76])
}
