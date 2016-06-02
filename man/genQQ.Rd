\name{genQQ}
\alias{genQQ}
\alias{generalizedQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalised quantile plot
}
\description{
Computes the empirical quantiles of the UH scores of a data vector and the theoretical quantiles of the standard exponential distribution. These quantiles are then plotted in a generalised QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
genQQ(data, gamma, plot = TRUE, main = "Generalised QQ-plot", ...)

generalizedQQ(data, gamma, plot = TRUE, main = "Generalised QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, typically Hill estimates are used.
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a generalised QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Generalised QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The \code{generalizedQQ} function is the same function but with a different name for compatibility with the old \code{S-Plus} code. 

The UH scores are defined as \eqn{UH_{j,n}=X_{n-j,n}H_{j,n}} with \eqn{H_{j,n}} the Hill estimates, but other positive estimates for the EVI can also be used. The appropriate positive estimates for the EVI need to be specified in \code{gamma}. The generalised QQ-plot then plots
\deqn{(\log((n+1)/(k+1)), \log(X_{n-k,n}H_{k,n}))}
for \eqn{k=1,\ldots,n-1}.

See Section 4.2.2 of Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{gqq.the}{Vector of the theoretical quantiles from a standard exponential distribution.}
\item{gqq.emp}{Vector of the empirical quantiles from the logarithm of the UH scores.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Beirlant, J., Vynckier, P. and Teugels, J.L. (1996). "Excess Function and Estimation of the Extreme-value Index." \emph{Bernoulli}, 2, 293--318.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ParetoQQ}}, \code{\link{Hill}}
}
\examples{
data(soa)

# Compute Hill estimator
H <- Hill(soa$size[1:5000], plot=FALSE)$gamma

# Generalised QQ-plot
genQQ(soa$size[1:5000], gamma=H)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
