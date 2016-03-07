\name{ParetoQQ}
\alias{ParetoQQ}
\title{
Pareto quantile plot
}
\description{
Computes the empirical quantiles of the log-transform of a data vector and the theoretical quantiles of the standard exponential distribution. These quantiles are then plotted in a Pareto QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
ParetoQQ(data, plot = TRUE, main = "Pareto QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a Pareto QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Pareto QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}

\details{
It can be easily seen that a log-transformed Pareto random variable is exponentially distributed.
We can hence obtain a Pareto QQ-plot from an exponential QQ-plot by replacing the empirical quantiles of the data vector by the empirical quantiles from the log-transformed data. We hence plot
\deqn{( -\log(1-i/(n+1)), \log X_{i,n} )}
for \eqn{i=1,...,n-1,}
with \eqn{X_{i,n}} the \eqn{i}-th order statistic of the data.

See Section 1.2 of Beirlant et al. (2004) for more details.
}

\value{
A list with following components:
\item{pqq.the}{Vector of the theoretical quantiles from a standard exponential distribution.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}

\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


\seealso{
\code{\link{ParetoQQ_der}}, \code{\link{ExpQQ}}, \code{\link{genQQ}}, \code{\link{LognormalQQ}}, \code{\link{WeibullQQ}}
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
