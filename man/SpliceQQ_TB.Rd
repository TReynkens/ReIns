\name{SpliceQQ_TB}
\alias{SpliceQQ_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing quantile plot using Turnbull estimator
}
\description{
This function plots the fitted quantile function of the spliced distribution versus quantiles based on the Turnbull survival function (which is suitable for right-censored or interval-censored data). 
}
\usage{
SpliceQQ_TB(L, U = L, censored, splicefit, 
            plot = TRUE, main = "Splicing QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
Vector of length \eqn{n} with the lower boundary of the intervals for interval censored data or the observed data for right censored data. In an insurance context, these are the paid amounts.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundary of the intervals. In an 
insurance context, these are the incurred values. By default, they are equal to \code{L}.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFiticPareto}}.
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
\details{
This QQ-plot is given by \deqn{(Q(p_j), \hat{Q}^{TB}]),} for \eqn{j=1,\ldots,n-1} where \eqn{Q} is the quantile function of the fitted splicing model, \eqn{\hat{Q}^{TB}} the empirical quantile function obtained using the Turnbull estimator and \eqn{p_j=j/(n+1)}. 

If the \pkg{interval} package is installed, the \code{\link[interval]{icfit}} function is used to compute the Turnbull estimator. Otherwise, \code{\link[survival]{survfit.formula}} from \pkg{survival} is used.

Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and right censored data should be entered as \code{L=trunclower} and \code{U=u}. The limits \code{trunclower} and \code{truncupper} are obtained from the \code{\link{SpliceFit}} object.

Use \code{\link{SpliceQQ}} for non-censored data. 
}
\value{
A list with following components:
\item{sqq.the}{Vector of the theoretical quantiles of the fitted spliced distribution.}
\item{sqq.emp}{Vector of the empirical quantiles from the data (based on the Turnbull estimator).}
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
\code{\link{SpliceQQ}}, \code{\link{qSplice}}, \code{\link{Turnbull}}, \code{\link[interval]{icfit}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceTB}}, \code{\link{SplicePP_TB}}, \code{\link{SpliceLL_TB}}
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
# Add 45 degree line
abline(0,1)
}
}

