\name{SplicePP}
\alias{SplicePP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PP-plot with fitted and empirical survival function
}
\description{
This function plots the fitted survival function of the spliced distribution versus the
empirical survival function (determined using the Empirical CDF (ECDF)). 
}
\usage{
SplicePP(X, splicefit, x = sort(X), log = FALSE, plot = TRUE, 
         main = "Splicing PP-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data used for fitting the distribution.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFitPareto}} or \code{\link{SpliceFitGPD}}.
}
  \item{x}{
Vector of points to plot the functions at. By default we plot them at the data points.
}
  \item{log}{
Logical indicating if minus the logarithms of the survival probabilities are plotted versus each other, default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the splicing PP-plot should be made, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Splicing PP-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The PP-plot consists of the points
\deqn{(1-\hat{F}(x_{i,n}), 1-\hat{F}_{spliced}(x_{i,n})))}
for \eqn{i=1,\ldots,n} with \eqn{n} the length of the data, \eqn{x_{i,n}} the \eqn{i}-th smallest observation,
\eqn{\hat{F}} the empirical distribution function and \eqn{\hat{F}_{spliced}} the fitted spliced distribution function.
The minus-log version of the PP-plot consists of
\deqn{(-\log(1-\hat{F}(x_{i,n})), -\log(1-\hat{F}_{spliced}(x_{i,n})))).}

Use \code{\link{SplicePP_TB}} for censored data. 

See Section 4.3 of Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{spp.the}{Vector of the theoretical probabilities \eqn{1-\hat{F}_{spliced}(x_{i,n})}.}
\item{spp.emp}{Vector of the empirical probabilities \eqn{1-\hat{F}(x_{i,n})}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}

\seealso{
\code{\link{SplicePP_TB}}, \code{\link{pSplice}}, \code{\link[stats]{ecdf}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitGPD}}, \code{\link{SpliceECDF}}, \code{\link{SpliceLL}}, \code{\link{SpliceQQ}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.6)



x <- seq(0, 20, 0.01)

# Plot of splicing CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of splicing PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")



# Fitted survival function and empirical survival function 
SpliceECDF(x, X, splicefit)

# Log-log plot with empirical survival function and fitted survival function
SpliceLL(x, X, splicefit)

# PP-plot of empirical survival function and fitted survival function
SplicePP(X, splicefit)

# PP-plot of empirical survival function and 
# fitted survival function with log-scales
SplicePP(X, splicefit, log=TRUE)

# Splicing QQ-plot
SpliceQQ(X, splicefit)
}
}
