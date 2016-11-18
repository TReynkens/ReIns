\name{cWeibullQQ}
\alias{cWeibullQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weibull quantile plot for right censored data
}
\description{
Weibull QQ-plot adapted for right censored data.
}
\usage{
cWeibullQQ(data, censored, plot = TRUE, main = "Weibull QQ-plot", ...)
}
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
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
The Weibull QQ-plot adapted for right censoring is given by
\deqn{( \log(-\log(1-F_{km}(Z_{j,n}))),  \log(Z_{j,n}) )}
for \eqn{j=1,\ldots,n-1,}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data and \eqn{F_{km}} the Kaplan-Meier estimator for the CDF.
Hence, it has the same empirical quantiles as an ordinary Weibull QQ-plot but replaces the theoretical quantiles \eqn{ \log(-\log(1-j/(n+1)))} by \eqn{\log(-\log(1-F_{km}(Z_{j,n})))}.

This QQ-plot is only suitable for right censored data.

In Beirlant et al. (2007), only a Pareto QQ-plot adapted for right-censored data is proposed. This QQ-plot is constructed using the same ideas, but is not described in the paper.
}
\value{
A list with following components:
\item{wqq.the}{Vector of the theoretical quantiles, see Details.}
\item{wqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Beirlant, J., Guillou, A., Dierckx, G. and Fils-Villetard, A. (2007). "Estimation of the Extreme Value Index and Extreme Quantiles Under Random Censoring." \emph{Extremes}, 10, 151--174.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{WeibullQQ}}, \code{\link{cExpQQ}}, \code{\link{cLognormalQQ}}, \code{\link{cParetoQQ}}, \code{\link{KaplanMeier}}
}
\examples{
# Set seed
set.seed(29072016)

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X, Y)

# Censoring indicator
censored <- (X>Y)

# Weibull QQ-plot adapted for right censoring
cWeibullQQ(Z, censored=censored)
}

