\name{ExcessEPD}
\alias{ExcessEPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using EPD estimates
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{M} and limit \eqn{L} using EPD estimates.
}
\usage{
ExcessEPD(data, gamma, delta, tau, M, L = Inf, warnings = TRUE, plot = TRUE, add = FALSE, 
          main = "Estimates for premium of excess-loss insurance", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, obtained from \code{\link{EPD}}.
}
  \item{delta}{
Vector of \eqn{n-1} estimates for \eqn{\delta}, obtained from \code{\link{EPD}}.
}
  \item{tau}{
Vector of \eqn{n-1} estimates for \eqn{\tau}, obtained from \code{\link{EPD}}.
}
  \item{M}{
The retention level of the (re-)insurance, default is \code{Inf}.
}
  \item{L}{
The limit of the (re-)insurance.
}
  \item{warnings}{
Logical indicating if warnings are displayed, default is \code{TRUE}.
}
 \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates for premium of excess-loss insurance"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
We need that \eqn{u \ge X_{n-k,n}}, the \eqn{(k+1)}-st largest observation.
If this is not the case, we return \code{NA} for the premium. A warning will be issued in
that case if \code{warnings=TRUE}.

The premium for the excess-loss insurance with retentation \eqn{M} and limit \eqn{L} is given by
\deqn{E(\min{(X-M)_+, L}) = I(M) - I(M+L)}
where \eqn{I(u)=\int_u^{\infty} (1-F(z)) dz} is the integrated tail function. When \eqn{L=\infty}, the premium is equal to \eqn{I(M)}.

We estimate the integrated tail function by
\deqn{ \hat{I}(u) = (k+1)/(n+1) * (X_{n-k,n})^{1/\hat{\gamma}} * ((1-\hat{\delta}/\hat{\gamma})(1/\hat{\gamma}-1)^{-1}u^{1/\hat{\gamma}} + 
\hat{\delta}/(\hat{\gamma}X_{n-k,n}^{\hat{\tau}})(1/\hat{\gamma}-\hat{\tau}-1)^{-1}u^{1+\hat{\tau}-1/\hat{\gamma}})}
with \eqn{\hat{\gamma}, \hat{\kappa}} and \eqn{\hat{\tau}} the estimates for the parameters of the EPD.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{premium}{The corresponding estimates for the premium.}
\item{M}{The retention level of the (re-)insurance.}
\item{L}{The limit of the (re-)insurance.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{EPD}}, \code{\link{ExcessHill}}, \code{\link{ExcessGPD}}
}
\examples{
data(secura)

# EPD estimator
epd <- EPD(secura$size)

# Premium of excess-loss insurance with retention M
M <- 10^7
ExcessEPD(secura$size,epd$gamma,epd$kappa,epd$tau,M=M,ylim=c(0,2*10^4))
}
