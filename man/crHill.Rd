\name{crHill}
\alias{crHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hill-type estimator for the conditional EVI
}
\description{
Hill-type estimator for the conditional Extreme Value Index (EVI) adapted for censored data.
}
\usage{
crHill(x, Xtilde, Ytilde, censored, h, 
       kernel = c("biweight", "normal", "uniform", "triangular", "epanechnikov"), 
       logk = FALSE, plot = FALSE, add = FALSE, main = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Value of the conditioning variable \eqn{X} to estimate the EVI at.
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
  \item{logk}{
Logical indicating if the Hill-type estimates are plotted as a function of \eqn{\log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{""} (no title).
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
This is a Hill-type estimator of the EVI of \eqn{Y} given \eqn{X=x}. The estimator uses the censored sample \eqn{(\tilde{X}_i, \tilde{Y}_i)}, for \eqn{i=1,\ldots,n}, where \eqn{X} and \eqn{Y} are censored at the same time. We assume that \eqn{Y} and the censoring variable are conditionally independent given \eqn{X}.

See Section 4.4.3 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding Hill-type estimates.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{crParetoQQ}}, \code{\link{crSurv}}, \code{\link{cHill}}
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

