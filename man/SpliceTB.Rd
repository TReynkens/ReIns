\name{SpliceTB}
\alias{SpliceTB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of fitted and Turnbull survival function
}
\description{
This function plots the fitted survival function of the splicing distribution together with the
Turnbull survival function (which is suitable for right-censored or interval-censored data). Moreover, \eqn{100\alpha\%} confidence bounds are added.
}
\usage{
SpliceTB(x, Z, I = Z, censored, splicefit, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to plot the functions at.
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
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitciPareto}}.
}
  \item{alpha}{
\eqn{100\alpha\%} is the confidence level for the confidence bounds. Default is \eqn{\alpha=0.05}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Use \code{\link{SpliceECDF}} for non-censored data. 
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm. \emph{Astin Bulletin}, 45, 729--758.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceECDF}}, \code{\link{pSplice}}, \code{\link{Turnbull}}, \code{\link{SpliceFitciPareto}}, \code{\link{SpliceLL_TB}}, \code{\link{SplicePP_TB}}
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
splicefit <- SpliceFitciPareto(Z=Z, censored=censored, const=0.9)

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
