\name{ExcessHill}
\alias{ExcessHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using a Pareto model
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{u} using a (truncated) Pareto model.
}
\usage{
ExcessHill(data, gamma, u, endpoint = Inf, warnings = TRUE, plot = TRUE, add = FALSE,
        main = "Estimates for premium of excess-loss insurance", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations.  
}
  \item{gamma}{
Vector of \eqn{n-1} estimates for the EVI, obtained from \code{\link{Hill}} or \code{\link{trHill}}.
}
  \item{u}{
The retention level of the (re-)insurance.
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
We need that \eqn{u \ge X_{n-k,n}}, the \eqn{(k+1)}-st largest observation.
If this is not the case, we return \code{NA} for the premium. A warning will be issued in
that case if \code{warnings=TRUE}. One should then use globale fits: \code{\link{ExcessSpliceHill}}.

We estimate the premium by
\deqn{ (k+1)/(n+1) / (1/H_{k,n}-1) * (X_{n-k,n})^{1/H_{k,n}} u^{1-H_{k,n}},}
with \eqn{H_{k,n}} the Hill estimator.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{premium}{The corresponding estimates for the premium.}
\item{u}{The retention level.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{ExcessEPD}}, \code{\link{ExcessGPD}}, \code{\link{ExcessSpliceHill}}
}
\examples{
data(secura)

# Hill estimator
H <- Hill(secura$size)

# Premium of excess-loss insurance with retention u
u <- 10^7
ExcessHill(secura$size,H$gamma,u=u)
}

