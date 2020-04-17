\name{ExpQQ}
\alias{ExpQQ}

\title{
Exponential quantile plot
}
\description{
Computes the empirical quantiles of a data vector and the theoretical quantiles of the standard exponential distribution. These quantiles are then plotted in an exponential QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
ExpQQ(data, plot = TRUE, main = "Exponential QQ-plot", ...)
}
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in an Exponential QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Exponential QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
The exponential QQ-plot is defined as
\deqn{( -\log(1-i/(n+1)), X_{i,n} )}
for \eqn{i=1,...,n,}
with \eqn{X_{i,n}} the \eqn{i}-th order statistic of the data.

Note that the mean excess plot is the derivative plot of the Exponential QQ-plot.

See Section 4.1 of Albrecher et al. (2017) for more details.
}

\value{
A list with following components:
\item{eqq.the}{Vector of the theoretical quantiles from a standard exponential distribution.}
\item{eqq.emp}{Vector of the empirical quantiles from the data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


\seealso{
\code{\link{MeanExcess}}, \code{\link{LognormalQQ}}, \code{\link{ParetoQQ}}, \code{\link{WeibullQQ}} 
}
\examples{
data(norwegianfire)

# Exponential QQ-plot for Norwegian Fire Insurance data for claims in 1976.
ExpQQ(norwegianfire$size[norwegianfire$year==76])

# Pareto QQ-plot for Norwegian Fire Insurance data for claims in 1976.
ParetoQQ(norwegianfire$size[norwegianfire$year==76])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

