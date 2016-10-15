\name{crParetoQQ}
\alias{crParetoQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional Pareto quantile plot for right censored data
}
\description{
Conditional Pareto QQ-plot adapted for right censored data.
}
\usage{
crParetoQQ(x, Xtilde, Ytilde, censored, h, 
           kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov"), 
           plot = TRUE, add = FALSE, main = "Pareto QQ-plot", type = "p", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Value of the conditioning variable \eqn{X} at which to make the conditional Pareto QQ-plot.
}
  \item{Xtilde}{
Vector of length \eqn{n} containing the censored sample of the conditioning variable \eqn{X}.
}
  \item{Ytilde}{
Vector of length \eqn{n} containing the censored sample  of the variable \eqn{Y}.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{h}{
Bandwidth of the non-parametric estimator for the conditional survival function (\code{\link{crSurv}}).
}
  \item{kernel}{
Kernel of the non-parametric estimator for the conditional survival function (\code{\link{crSurv}}). One of "biweight" (default), "normal", "uniform", "triangular" and "epanechnikov".
}
 \item{plot}{
Logical indicating if the quantiles should be plotted in a Pareto QQ-plot, default is \code{TRUE}.
}
  \item{add}{
Logical indicating if the quantiles should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Pareto QQ-plot"}.
}
  \item{type}{
Type of the plot, default is \code{"p"} meaning points are plotted, see \code{\link[graphics]{plot}} for more details.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
We construct a Pareto QQ-plot for \eqn{Y} conditional on \eqn{X=x} using the censored sample \eqn{(\tilde{X}_i, \tilde{Y}_i)}, for \eqn{i=1,\ldots,n}, where \eqn{X} and \eqn{Y} are censored at the same time. We assume that \eqn{Y} and the censoring variable are conditionally independent given \eqn{X}.

The conditional Pareto QQ-plot adapted for right censoring is given by
\deqn{( -\log(1-\hat{F}_{Y|X}(\tilde{Y}_{j,n}|x)), \log \tilde{Y}_{j,n} )}
for \eqn{j=1,\ldots,n-1,}
with \eqn{\tilde{Y}_{i,n}} the \eqn{i}-th order statistic of the censored data and \eqn{\hat{F}_{Y|X}(y|x)} the non-parametric estimator for the conditional CDF of Akritas and Van Keilegom (2003), see \code{\link{crSurv}}.

See Section 4.4.3 in Albrecher et al. (2017) for more details.
}

\value{
A list with following components:
\item{pqq.the}{Vector of the theoretical quantiles, see Details.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed \eqn{Y} data.}
}
\references{
Akritas, M.G. and Van Keilegom, I. (2003). "Estimation of Bivariate and Marginal Distributions
With Censored Data." \emph{Journal of the Royal Statistical Society: Series B}, 65, 457--471.

Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link{crSurv}}, \code{\link{crHill}}, \code{\link{cParetoQQ}}
}
\examples{
# Set seed
set.seed(29072016)

# Pareto random sample
Y <- rpareto(200, shape=2)

# Censoring variable
C <- rpareto(200, shape=1)

# Observed (censored) sample of variable Y
Ytilde <- pmin(Y, C)

# Censoring indicator
censored <- (Y>C)

# Conditioning variable
X <- seq(1, 10, length.out=length(Y))

# Observed (censored) sample of conditioning variable
Xtilde <- X
Xtilde[censored] <- X[censored] - runif(sum(censored), 0, 1)


# Conditional Pareto QQ-plot
crParetoQQ(x=1, Xtilde=Xtilde, Ytilde=Ytilde, censored=censored, h=2)

# Plot Hill-type estimates
crHill(x=1, Xtilde, Ytilde, censored, h=2, plot=TRUE)
}

