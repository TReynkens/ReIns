\name{SplicePP_TB}
\alias{SplicePP_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PP-plot with fitted and Turnbull survival function
}
\description{
This function plots the fitted survival function of the splicing distribution versus the Turnbull survival function (which is suitable for right-censored or interval-censored data). 
}
\usage{
SplicePP_TB(x = sort(Z), Z, I = Z, censored, splicefit, log = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
Vector of points to plot the functions at. By default we plot it at the data points \code{Z}.
}
  \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals for interval censored data or the observed data for right censored data. In an insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values. By default, they are equal to \code{Z} which means that we consider right-censored data. 
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFitcHill}}.
}
  \item{log}{
Logical indicating if the logarithms of the survival probabilities are plotted versus each other, default is \code{FALSE}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Use \code{\link{SplicePP}} for non-censored data. 
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SplicePP}}, \code{\link{pSplice}}, \code{\link{Turnbull}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceTB}}, \code{\link{SpliceLL_TB}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y  <- rpareto(500, shape=1)

# Observed sample
Z <- pmin(X,Y)

# Censoring indicator
censored <- (X>Y)

# Splice ME and Pareto
splicefit <- SpliceFitcHill(Z=Z, censored=censored, const=0.9)



x <- seq(0,10^2,0.1)

# Plot of splicing CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of splicing PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")


# Fitted survival function and Turnbull survival function 
SpliceTB(x, Z, censored=censored, splicefit=splicefit)

# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB(x, Z, censored=censored, splicefit=splicefit)

# PP-plot of Turnbull survival function and fitted survival function
SplicePP_TB(x, Z, censored=censored, splicefit=splicefit)

# PP-plot of Turnbull survival function and 
# fitted survival function with log-scales
SplicePP_TB(x, Z, censored=censored, splicefit=splicefit, log=TRUE)

}
}
