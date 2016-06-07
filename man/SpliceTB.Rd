\name{SpliceTB}
\alias{SpliceTB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of fitted and Turnbull survival function
}
\description{
This function plots the fitted survival function of the spliced distribution together with the
Turnbull survival function (which is suitable for interval censored data). Moreover, \eqn{100\alpha\%} confidence bounds are added.
}
\usage{
SpliceTB(x = sort(L), L, U = L, censored, splicefit, alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to plot the functions at. By default we plot it at the points \code{L}.
}
  \item{L}{
Vector of length \eqn{n} with the lower boundary of the intervals for interval censored data or the observed data for right censored data.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundary of the intervals. By default, they are equal to \code{L}.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFiticPareto}}.
}
  \item{alpha}{
\eqn{100\alpha\%} is the confidence level for the confidence bounds. Default is \eqn{\alpha=0.05}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and left censored data should be entered as \code{L=trunclower} and \code{U=u}. The limits \code{trunclower} and \code{truncupper} are obtained from the \code{\link{SpliceFit}} object.

Use \code{\link{SpliceECDF}} for non-censored data. 

See Section 4.3 of Albrecher et al. (2016) for more details.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceECDF}}, \code{\link{pSplice}}, \code{\link{Turnbull}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceLL_TB}}, \code{\link{SplicePP_TB}}, \code{\link{SpliceQQ_TB}}
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

# Right boundary
U <- Z
U[censored] <- Inf

# Splice ME and Pareto
splicefit <- SpliceFiticPareto(L=Z, U=U, censored=censored, tsplice=quantile(Z,0.9))



x <- seq(0,20,0.1)

# Plot of splicing CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of splicing PDF
plot(x, dSplice(x, splicefit), type="l", xlab="x", ylab="f(x)")


# Fitted survival function and Turnbull survival function 
SpliceTB(x, L=Z, U=U, censored=censored, splicefit=splicefit)


# Log-log plot with Turnbull survival function and fitted survival function
SpliceLL_TB(x, L=Z, U=U, censored=censored, splicefit=splicefit)


# PP-plot of Turnbull survival function and fitted survival function
SplicePP_TB(L=Z, U=U, censored=censored, splicefit=splicefit)

# PP-plot of Turnbull survival function and 
# fitted survival function with log-scales
SplicePP_TB(L=Z, U=U, censored=censored, splicefit=splicefit, log=TRUE)

# QQ-plot using Turnbull survival function and fitted survival function
SpliceQQ_TB(L=Z, U=U, censored=censored, splicefit=splicefit)
}
}
