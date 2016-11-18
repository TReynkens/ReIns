\name{crSurv}
\alias{crSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-parametric estimator of conditional survival function
}
\description{
Non-parametric estimator of the conditional survival function of \eqn{Y} given \eqn{X} for censored data, see Akritas and Van Keilegom (2003).
}
\usage{
crSurv(x, y, Xtilde, Ytilde, censored, h, 
       kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The value of the conditioning variable \eqn{X} to evaluate the survival function at. \code{x} needs to be a single number or a vector with the same length as \code{y}.
}
  \item{y}{
The value(s) of the variable \eqn{Y} to evaluate the survival function at.
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
Bandwidth of the non-parametric estimator.
}
  \item{kernel}{
Kernel of the non-parametric estimator. One of \code{"biweight"} (default), \code{"normal"}, \code{"uniform"}, \code{"triangular"} and \code{"epanechnikov"}.
}
}
\details{
We estimate the conditional survival function \deqn{1-F_{Y|X}(y|x)} using the censored sample \eqn{(\tilde{X}_i, \tilde{Y}_i)}, for \eqn{i=1,\ldots,n}, where \eqn{X} and \eqn{Y} are censored at the same time. We assume that \eqn{Y} and the censoring variable are conditionally independent given \eqn{X}.

The estimator is given by
\deqn{1-\hat{F}_{Y|X}(y|x) = \prod_{\tilde{Y}_i \le y} (1-W_{n,i}(x;h_n)/(\sum_{j=1}^nW_{n,j}(x;h_n) I\{\tilde{Y}_j \ge \tilde{Y}_i\}))^{\Delta_i}}
where \eqn{\Delta_i=1} when \eqn{(\tilde{X}_i, \tilde{Y}_i)} is censored and 0 otherwise. The weights are given by
\deqn{W_{n,i}(x;h_n) = K((x-\tilde{X}_i)/h_n)/\sum_{\Delta_j=1}K((x-\tilde{X}_j)/h_n)} when \eqn{\Delta_i=1} and 0 otherwise.

See Section 4.4.3 in Albrecher et al. (2017) for more details.
}
\value{
Estimates for \eqn{1-F_{Y|X}(y|x)} as described above.
}
\references{
Akritas, M.G. and Van Keilegom, I. (2003). "Estimation of Bivariate and Marginal Distributions
With Censored Data." \emph{Journal of the Royal Statistical Society: Series B}, 65, 457--471.

Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{crParetoQQ}}, \code{\link{crHill}}
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

# Plot estimates of the conditional survival function
x <- 5
y <- seq(0, 5, 1/100)
plot(y, crSurv(x, y, Xtilde=Xtilde, Ytilde=Ytilde, censored=censored, h=5), type="l", 
     xlab="y", ylab="Conditional survival function")
}