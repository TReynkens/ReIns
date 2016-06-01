\name{WeibullQQ}
\alias{WeibullQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weibull quantile plot
}
\description{
Computes the empirical quantiles of the log-transform of a data vector and the theoretical quantiles of the standard Weibull distribution. These quantiles are then plotted in a Weibull QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
WeibullQQ(data, plot = TRUE, main = "Weibull QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a Weibull QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Weibull QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The Weibull QQ-plot is given by
\deqn{(\log(-\log(1-i/(n+1))),\log X_{i,n})}
for \eqn{i=1,...,n,} with \eqn{X_{i,n}} the \eqn{i}-th order statistic of the data.

See Section 4.1 of Albrecher et al. (2016) for more details.
}

\value{
A list with following components:
\item{wqq.the}{Vector of the theoretical quantiles from a standard Weibull distribution.}
\item{wqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}

\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{WeibullQQ_der}}, \code{\link{ExpQQ}}, \code{\link{LognormalQQ}}, \code{\link{ParetoQQ}}
}
\examples{
data(norwegianfire)

# Weibull QQ-plot for Norwegian Fire Insurance data for claims in 1976.
WeibullQQ(norwegianfire$size[norwegianfire$year==76])

# Derivative of Weibull QQ-plot for Norwegian Fire Insurance data for claims in 1976.
WeibullQQ_der(norwegianfire$size[norwegianfire$year==76])
}