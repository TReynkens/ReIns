\name{ExcessPareto}
\alias{ExcessPareto}
\alias{ExcessHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using a Pareto model
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{R} and limit \eqn{L} using a (truncated) Pareto model.
}
\usage{
ExcessPareto(data, gamma, R, L = Inf, endpoint = Inf, warnings = TRUE, plot = TRUE, 
        add = FALSE, main = "Estimates for premium of excess-loss insurance", ...)
        
ExcessHill(data, gamma, R, L = Inf, endpoint = Inf, warnings = TRUE, plot = TRUE, 
        add = FALSE, main = "Estimates for premium of excess-loss insurance", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, obtained from \code{\link{Hill}} or \code{\link{trHill}}.
}
 \item{R}{
The retention level of the (re-)insurance.
}
  \item{L}{
The limit of the (re-)insurance, default is \code{Inf}.
}
\item{endpoint}{
Endpoint for the truncated Pareto distribution. When \code{Inf}, the default, the ordinary Pareto model is used.
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
We need that \eqn{u \ge X_{n-k,n}}, the \eqn{(k+1)}-th largest observation.
If this is not the case, we return \code{NA} for the premium. A warning will be issued in
that case if \code{warnings=TRUE}. One should then use global fits: \code{\link{ExcessSplice}}.

The premium for the excess-loss insurance with retention \eqn{R} and limit \eqn{L} is given by
\deqn{E(\min{(X-R)_+, L}) = \Pi(R) - \Pi(R+L)}
where \eqn{\Pi(u)=E((X-u)_+)=\int_u^{\infty} (1-F(z)) dz} is the premium of the excess-loss insurance with retention \eqn{u}. When \eqn{L=\infty}, the premium is equal to \eqn{\Pi(R)}.

We estimate \eqn{\Pi} (for the untruncated Pareto distribution) by
\deqn{ \hat{\Pi}(u) = (k+1)/(n+1) / (1/H_{k,n}-1) * (X_{n-k,n}^{1/H_{k,n}} u^{1-1/H_{k,n}}),}
with \eqn{H_{k,n}} the Hill estimator.

The \code{ExcessHill} function is the same function but with a different name for compatibility with old versions of the package.

See Section 4.6 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{premium}{The corresponding estimates for the premium.}
\item{R}{The retention level of the (re-)insurance.}
\item{L}{The limit of the (re-)insurance.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{ExcessEPD}}, \code{\link{ExcessGPD}}, \code{\link{ExcessSplice}}
}
\examples{
data(secura)

# Hill estimator
H <- Hill(secura$size)

# Premium of excess-loss insurance with retention R
R <- 10^7
ExcessHill(secura$size, H$gamma, R=R)
}

