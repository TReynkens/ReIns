\name{trMLE}
\alias{trMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MLE estimator for upper truncated data
}
\description{
Computes the ML estimator for the extreme value index, adapted for upper truncation, as a function of the tail parameter \eqn{k} (Beirlant et al., 2017). Optionally, these estimates are plotted as a function of \eqn{k}.
}
\usage{
trMLE(data, start = c(1, 1), eps = 10^(-10), 
      plot = TRUE, add = FALSE, main = "Estimates for EVI", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 Vector of \eqn{n} observations. 
}
  \item{start}{
  Starting values for \eqn{\gamma} and \eqn{\tau} for the numerical optimisation.
}
  \item{eps}{
Numerical tolerance, see Details. By default it is equal to \code{10^(-10)}.
}
 \item{plot}{
Logical indicating if the estimates of \eqn{\gamma} should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates of \eqn{\gamma} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of the EVI"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
We compute the MLE for the \eqn{\gamma} and \eqn{\sigma} parameters of the truncated GPD.
For numerical reasons, we compute the MLE for \eqn{\tau=\gamma/\sigma} and transform this estimate to \eqn{\sigma}.

The log-likelihood is given by
\deqn{(k-1) \ln \tau - (k-1) \ln \xi- ( 1 + 1/\xi)\sum_{j=2}^k \ln (1+\tau E_{j,k}) -(k-1) \ln( 1- (1+ \tau E_{1,k})^{-1/\xi})}
with \eqn{E_{j,k} = X_{n-j+1,n}-X_{n-k,n}}.

In order to meet the restrictions \eqn{\sigma=\xi/\tau>0} and  \eqn{1+\tau E_{j,k}>0} for \eqn{j=1,\ldots,k}, we require the estimates of these quantities to be larger than the numerical tolerance value \code{eps}.

See Beirlant et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{gamma}{Vector of the corresponding estimates for \eqn{\gamma}.}
\item{tau}{Vector of the corresponding estimates for \eqn{\tau}.}
\item{sigma}{Vector of the corresponding estimates for \eqn{\sigma}.}
\item{conv}{Convergence indicator of \code{\link[stats]{optim}}.}
}
\references{
Beirlant, J., Fraga Alves, M. I. and Reynkens, T. (2017). "Fitting Tails Affected by Truncation". \emph{Electronic Journal of Statistics}, 11(1), 2026--2065.
}
\author{
Tom Reynkens.
}
\seealso{
\code{\link{trDTMLE}}, \code{\link{trEndpointMLE}}, \code{\link{trProbMLE}}, \code{\link{trQuantMLE}}, \code{\link{trTestMLE}}, \code{\link{trHill}}, \code{\link{GPDmle}}
}
\examples{
# Sample from GPD truncated at 99\% quantile
gamma <- 0.5
sigma <- 1.5
X <- rtgpd(n=250, gamma=gamma, sigma=sigma, endpoint=qgpd(0.99, gamma=gamma, sigma=sigma))

# Truncated ML estimator
trmle <- trMLE(X, plot=TRUE, ylim=c(0,2))
}

