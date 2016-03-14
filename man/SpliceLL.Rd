\name{SpliceLL}
\alias{SpliceLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
LL-plot with fitted and empirical survival function
}
\description{
This function plots the logarithm of the empirical survival function (determined using the Empirical CDF (ECDF)) versus the logarithm of the data.
Moreover, the logarithm of the fitted survival function of the spliced distribution is added.
}
\usage{
SpliceLL(x = sort(X), X, splicefit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to plot the fitted survival function at. By default we plot it at the data points.
}
  \item{X}{
Data used for fitting the distribution.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFitPareto}} or \code{\link{SpliceFitGPD}}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Use \code{\link{SpliceLL_TB}} for censored data. 
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
\code{\link{SpliceLL_TB}}, \code{\link{pSplice}}, \code{\link[stats]{ecdf}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitGPD}}, \code{\link{SpliceECDF}}, \code{\link{SplicePP}}, \code{\link{SpliceQQ}}
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
