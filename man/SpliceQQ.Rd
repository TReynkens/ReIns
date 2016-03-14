\name{SpliceQQ}
\alias{SpliceQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing quantile plot
}
\description{
Computes the empirical quantiles of a data vector and the theoretical quantiles of the fitted spliced distribution. These quantiles are then plotted in a splicing QQ-plot with the theoretical quantiles on the \eqn{x}-axis and the empirical quantiles on the \eqn{y}-axis.
}
\usage{
SpliceQQ(X, splicefit, p = NULL, plot = TRUE, main = "Splicing QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  Vector of \eqn{n} observations. 
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFitPareto}} or \code{\link{SpliceFitGPD}}.
}
 \item{p}{
  Vector of probabilities used in the QQ-plot. If \code{NULL}, the default, we take \code{p} equal to \code{1/(n+1),...,n/(n+1)}.
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a splicing QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Splicing QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\value{
A list with following components:
\item{sqq.the}{Vector of the theoretical quantiles of the fitted spliced distribution.}
\item{sqq.emp}{Vector of the empirical quantiles from the data.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm. \emph{Astin Bulletin}, 45, 729--758.
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{SpliceQQ_TB}}, \code{\link{qSplice}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitGPD}}, \code{\link{SpliceECDF}}, \code{\link{SpliceLL}}, \code{\link{SplicePP}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)



x <- seq(0,20,0.01)

# Plot of splicing CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of splicing PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")



# Fitted survival function and empirical survival function 
SpliceECDF(x, X, splicefit)

# Log-log plot with empirical survival function and fitted survival function
SpliceLL(x, X, splicefit)

# PP-plot of empirical survival function and fitted survival function
SplicePP(x, X, splicefit)

# PP-plot of empirical survival function and 
# fitted survival function with log-scales
SplicePP(x, X, splicefit, log=TRUE)

# Splicing QQ-plot
SpliceQQ(X, splicefit)

}
}
