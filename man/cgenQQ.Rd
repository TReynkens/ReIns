\name{cgenQQ}
\alias{cgenQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalised QQplot for right censored data
}
\description{
Generalised QQplot adapted for right censored data.
}
\usage{
cgenQQ(data, censored, plot = TRUE, main = "Generalised QQ-plot", ...)
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
Logical indicating if the quantiles should be plotted in a generalised QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Generalised QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The generalised QQplot adapted for right censoring is given by
\deqn{( -\log(1-F_{km}(Z_{n-j+1,n})), \log UH_{n-j+1,n} )}
for \eqn{j=1,...,n-1,} with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{UH_{n-j+1,n} = Z_{n-j+1,n} H_{j-1,n}},  \eqn{H_{j-1,n}^c} the (ordinary) Hill estimator and \eqn{F_{km}} the Kaplan-Meier estimator for the CDF.

This estimator is only suitable for right censored data!
}
\value{
A list with following components:
\item{gqq.the}{Vector of the theoretical quantiles, see Details.}
\item{gqq.emp}{Vector of the empirical quantiles from the UH scores, see Details.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{genQQ}}, \code{\link{cgenHill}}, \code{\link{cParetoQQ}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Generalised QQplot adapted for right censoring
cgenQQ(secura$size,censored=censored)
}