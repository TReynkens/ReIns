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
SpliceLL(x = sort(X), X, splicefit, plot = TRUE, main = "Splicing LL-plot", ...)
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
  \item{plot}{
Logical indicating if the splicing LL-plot should be made, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"Splicing LL-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The LL-plot consists of the points
\deqn{(\log(x_{i,n}), \log(1-\hat{F}(x_{i,n})))}
for \eqn{i=1,\ldots,n} with \eqn{n} the length of the data, \eqn{x_{i,n}} the \eqn{i}-th smallest observation
and \eqn{\hat{F}} the empirical distribution function.
Then, the line 
\deqn{(\log(x), \log(1-\hat{F}_{spliced}(x))),}
with \eqn{\hat{F}_{spliced}} the fitted spliced distribution function, is added.

Use \code{\link{SpliceLL_TB}} for censored data.

See Reynkens et al. (2017) and Section 4.3.1 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{logX}{Vector of the logarithms of the sorted data.}
\item{sll.the}{Vector of the theoretical log-probabilities \eqn{\log(1-\hat{F}_{spliced}(x))}.}
\item{logx}{Vector of the logarithms of the points to plot the fitted survival function at.}
\item{sll.emp}{Vector of the empirical log-probabilities \eqn{\log(1-\hat{F}(x_{i,n}))}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{Insurance: Mathematics and Economics}, 77, 65--77.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
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



x <- seq(0, 20, 0.01)

# Plot of spliced CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of spliced PDF
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
