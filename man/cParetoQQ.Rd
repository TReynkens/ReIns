\name{cParetoQQ}
\alias{cParetoQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pareto quantile plot for right censored data
}
\description{
Pareto QQ-plot adapted for right censored data.
}
\usage{
cParetoQQ(data, censored, plot = TRUE, main = "Pareto QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
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
The Pareto QQ-plot adapted for right censoring is given by
\deqn{( -\log(1-F_{km}(Z_{j,n})), \log Z_{j,n} )}
for \eqn{j=1,\ldots,n-1,}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data and \eqn{F_{km}} the Kaplan-Meier estimator for the CDF.
Hence, it has the same empirical quantiles as an ordinary Pareto QQ-plot but replaces the theoretical quantiles \eqn{-\log(1-j/(n+1))} by \eqn{-\log(1-F_{km}(Z_{j,n}))}.

This QQ-plot is only suitable for right censored data, use \code{link{ciParetQQ}} for interval censored data.
}

\value{
A list with following components:
\item{pqq.the}{Vector of the theoretical quantiles, see Details.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Beirlant, J., Guillou, A., Dierckx, G. and Fils-Villetard, A. (2007). "Estimation of the Extreme Value Index and Extreme Quantiles Under Random Censoring." \emph{Extremes}, 10, 151--174.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ParetoQQ}}, \code{cExpQQ}, \code{cLognormalQQ}, \code{\link{cWeibullQQ}}, \code{\link{cHill}}, \code{\link{KaplanMeier}}
}
\examples{
# Set seed
set.seed(29072016)

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y  <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X, Y)

# Censoring indicator
censored <- (X>Y)

# Pareto QQ-plot adapted for right censoring
cParetoQQ(Z, censored=censored)
}
