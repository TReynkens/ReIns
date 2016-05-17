\name{WeibullQQ_der}
\alias{WeibullQQ_der}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Derivative plot of the Weibull QQ-plot
}
\description{
Computes the derivative plot of the Weibull QQ-plot. These values can be plotted as a function of the data or as a function of the tail parameter \eqn{k}.
}
\usage{
WeibullQQ_der(data, k = FALSE, plot = TRUE, 
              main = "Derivative plot of Weibull QQ-plot", ...)
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
Title for the plot, default is \code{"Derivative plot of Weibull QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The derivative plot of a Weibull QQ-plot is \deqn{(k, H_{k,n}/W_{k,n})} or \deqn{(\log X_{n-k,n}, H_{k,n}/W_{k,n})}
with \eqn{H_{k,n}} the Hill estimates and
\deqn{W_{k,n}=1/k*\sum_{j=1}^k \log(\log((n+1)/j))  -  \log(\log((n+1)/(k+1))).}
}

\value{
A list with following components:
\item{xval}{Vector of the x-values of the plot (\eqn{k} or \eqn{\log X_{n-k,n}}).}
\item{yval}{Vector of the derivative values.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{WeibullQQ}}, \code{\link{MeanExcess}}, \code{\link{LognormalQQ_der}}, \code{\link{ParetoQQ_der}}
}
\examples{
data(norwegianfire)

# Weibull QQ-plot for Norwegian Fire Insurance data for claims in 1976.
WeibullQQ(norwegianfire$size[norwegianfire$year==76])

# Derivative of Weibull QQ-plot for Norwegian Fire Insurance data for claims in 1976.
WeibullQQ_der(norwegianfire$size[norwegianfire$year==76])
}