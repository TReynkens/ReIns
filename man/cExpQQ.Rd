\name{cExpQQ}
\alias{cExpQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exponential QQ-plot for right censored data
}
\description{
Exponential QQ-plot adapted for right censored data.
}
\usage{
cExpQQ(data, censored, plot = TRUE, main = "Exponential QQ-plot", ...)
}
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in an exponential QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Exponential QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The exponential QQ-plot adapted for right censoring is given by
\deqn{( -\log(1-F_{km}(Z_{n-j+1,n})),  Z_{n-j+1,n} )}
for \eqn{j=1,...,n-1,}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data and \eqn{F_{km}} the Kaplan-Meier estimator for the CDF.
Hence, it has the same empirical quantiles as an ordinary exponential QQ-plot but replaces the theoretical quantiles \eqn{-\log(j/(n+1))} by \eqn{-\log(1-F_{km}(Z_{n-j+1,n}))}.

This estimator is only suitable for right censored data.
}

\value{
A list with following components:
\item{eqq.the}{Vector of the theoretical quantiles, see Details.}
\item{eqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ExpQQ}}, \code{\link{cLognormalQQ}}, \code{\link{cParetoQQ}}, \code{\link{KaplanMeier}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Exponential QQ-plot adapted for right censoring
cExpQQ(secura$size,censored=censored)
}
