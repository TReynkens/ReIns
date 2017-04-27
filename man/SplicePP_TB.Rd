\name{SplicePP_TB}
\alias{SplicePP_TB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PP-plot with fitted and Turnbull survival function
}
\description{
This function plots the fitted survival function of the spliced distribution versus the Turnbull survival function (which is suitable for interval censored data). 
}
\usage{
SplicePP_TB(L, U = L, censored, splicefit, x = NULL, log = FALSE, plot = TRUE,
            main = "Splicing PP-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
Vector of length \eqn{n} with the lower boundaries of the intervals for interval censored data or the observed data for right censored data.
}
  \item{U}{
Vector of length \eqn{n} with the upper boundaries of the intervals. By default, they are equal to \code{L}.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object, e.g. output from \code{\link{SpliceFiticPareto}}.
}
 \item{x}{
Vector of points to plot the functions at. When \code{NULL}, the default, the empirical quantiles for \eqn{1/(n+1), \ldots, n/(n+1)}, obtained using the Turnbull estimator, are used.
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
\deqn{(1-\hat{F}^{TB}(x_i), 1-\hat{F}_{spliced}(x_i)))}
for \eqn{i=1,\ldots,n} with \eqn{n} the length of the data, \eqn{x_i=\hat{Q}^{TB}(p_i)} where \eqn{p_i=i/(n+1)}, \eqn{\hat{Q}^{TB}} is the quantile function obtained using the Turnbull estimator,
\eqn{\hat{F}^{TB}} the Turnbull estimator for the distribution function and \eqn{\hat{F}_{spliced}} the fitted spliced distribution function.
The minus-log version of the PP-plot consists of
\deqn{(-\log(1-\hat{F}^{TB}(x_i)), -\log(1-\hat{F}_{spliced}(x_i))).}



Right censored data should be entered as \code{L=l} and \code{U=truncupper}, and left censored data should be entered as \code{L=trunclower} and \code{U=u}. The limits \code{trunclower} and \code{truncupper} are obtained from the \code{\link{SpliceFit}} object.

If the \pkg{interval} package is installed, the \code{\link[interval]{icfit}} function is used to compute the Turnbull estimator. Otherwise, \code{\link[survival]{survfit.formula}} from \pkg{survival} is used.

Use \code{\link{SplicePP}} for non-censored data. 

See Reynkens et al. (2017) and Section 4.3.2 in Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{spp.the}{Vector of the theoretical probabilities \eqn{1-\hat{F}_{spliced}(x_i)} (when \code{log=FALSE}) or \eqn{-\log(1-\hat{F}_{spliced}(x_i))} (when \code{log=TRUE}).}
\item{spp.emp}{Vector of the empirical probabilities \eqn{1-\hat{F}^{TB}(x_i)} (when \code{log=FALSE}) or \eqn{-\log(1-\hat{F}^{TB}(x_i))} (when \code{log=TRUE}).}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SplicePP}}, \code{\link{pSplice}}, \code{\link{Turnbull}}, \code{\link[interval]{icfit}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceTB}}, \code{\link{SpliceLL_TB}}, \code{\link{SpliceQQ_TB}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(500, shape=2)

# Censoring variable
Y <- rpareto(500, shape=1)

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

# Plot of spliced CDF
plot(x, pSplice(x, splicefit), type="l", xlab="x", ylab="F(x)")

# Plot of spliced PDF
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
