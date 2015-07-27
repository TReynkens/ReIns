\name{ciParetoQQ}
\alias{ciParetoQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pareto QQ-plot for interval censored data
}
\description{
Pareto QQ-plot adapted for interval censored data.
}
\usage{
ciParetoQQ(Z, I, censored, plot = TRUE, main = "Pareto QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals. In an 
insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values.
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
The Pareto QQ-plot adapted for interval censoring is given by
\deqn{( -\log(1-F_{TB}(Z_{n-j+1,n})), \log Z_{n-j+1,n} )}
for \eqn{j=1,...,n-1,}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data and \eqn{F_{TB}} the Turnbull estimator for the CDF.
Hence, it has the same empirical quantiles as an ordinary Pareto QQ-plot but replaces the theoretical quantiles \eqn{-\log(j/(n+1))} by \eqn{-\log(1-F_{TB}(Z_{n-j+1,n}))}.
}
\value{
A list with following components:
\item{pqq.the}{Vector of the theoretical quantiles, see Details.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ParetoQQ}}, \code{\link{cParetoQQ}}, \code{\link{ciExpQQ}}, \code{ciLognormalQQ}, \code{\link{ciHill}}, \code{\link{Turnbull}}
}
\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Incurred values
I <- secura$size+rchisq(length(secura$size),2)

# Pareto QQ-plot adapted for interval censoring
ciParetoQQ(secura$size,I,censored=censored)
}
