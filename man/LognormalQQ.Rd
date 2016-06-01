\name{LognormalQQ}
\alias{LognormalQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-normal quantile plot 
}
\description{
Computes the empirical quantiles of the log-transform of a data vector and the theoretical quantiles of the standard normal distribution. These quantiles are then plotted in a log-normal QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
LognormalQQ(data, plot = TRUE, main = "Log-normal QQ-plot", ...)
}
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a log-normal QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Log-normal QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
By definition, a log-transformed log-normal random variable is normally distributed.
We can thus obtain a log-normal QQ-plot from a normal QQ-plot by replacing the empirical quantiles of the data vector by the empirical quantiles from the log-transformed data. We hence plot
\deqn{(\Phi^{-1}(i/(n+1)), \log(X_{i,n}) )} for \eqn{i=1,\ldots,n,} where \eqn{\phi} is the normal PDF and \eqn{\Phi} is the normal CDF.

See Section 4.1 of Albrecher et al. (2016) for more details.
}

\value{
A list with following components:
\item{lnqq.the}{Vector of the theoretical quantiles from a standard normal distribution.}
\item{lnqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ExpQQ}}, \code{\link{ParetoQQ}}, \code{\link{WeibullQQ}} 
}
\examples{
data(norwegianfire)

# Log-normal QQ-plot for Norwegian Fire Insurance data for claims in 1976.
LognormalQQ(norwegianfire$size[norwegianfire$year==76])


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


