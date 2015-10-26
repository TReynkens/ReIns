\name{SpliceECDF}
\alias{SpliceECDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of fitted and empirical survival function
}
\description{
This function plots the fitted survival function of the splicing distribution together with the
empirical survival function (determined using the Empirical CDF (ECDF)). Moreover, \eqn{100\alpha\%} confidence bounds are added.
}
\usage{
SpliceECDF(x, X, splicefit, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to plot the functions at.
}
  \item{X}{
Data used for fitting the distribution.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFitPareto}} or \code{\link{SpliceFitGPD}}.
}
  \item{alpha}{
\eqn{100\alpha\%} is the confidence level for the confidence bounds.  Default is \eqn{\alpha=0.05}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Use \code{\link{SpliceTB}} for censored data. 
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm. \emph{Astin Bulletin}, 45, 729--758.
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{SpliceTB}}, \code{\link{pSplice}}, \code{\link[stats]{ecdf}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitGPD}}, \code{\link{SpliceLL}}, \code{\link{SplicePP}}, \code{\link{SpliceQQ}}
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
