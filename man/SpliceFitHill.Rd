\name{SpliceFitHill}
\alias{SpliceFitHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing of mixed Erlang and Pareto
}
\description{
Fit spliced distribution of a mixed Erlang and Pareto distribution(s). The shape parameter of the Pareto distribution(s) is determined using the Hill estimator. 
}
\usage{
SpliceFitHill(X, const, M = 3, s = 1:10, trunclower = 0, 
              EVTtruncation = FALSE, ncores = NULL, criterium = c("BIC","AIC"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data used for fitting the distribution.
}
  \item{const}{
Vector of length \eqn{l} containing the probabilities of the quantiles where the distributions will be spliced (splicing points). The ME distribution will be spliced with \eqn{l} Pareto distributions.
}

  \item{M}{
Initial number of Erlang mixtures, default is 3. This number can change when determining
an optimal mixed Erlang fit using an AIC criterion.
}
  \item{s}{
Vector of spread factors for the EM algorithm, default is \code{1:10}. We loop over these factors when determining
an optimal mixed Erlang fit using an AIC criterion.
}
  \item{trunclower}{
Lower truncation point. Default is 0.
}
  \item{EVTtruncation}{
Logical vector of length \eqn{l} indicating for each Pareto distribution if it is a truncated Pareto distribution.
A single logical can also be provided which will then be used for all \eqn{l} distributions.
}
  \item{ncores}{
Number of cores to use when determining an optimal mixed Erlang fit using an AIC criterion.
When \code{NULL}, \code{max(nc-1,1)} cores are used where \code{nc} is the number of cores as determined by \code{\link{detectCores}}.
}
\item{criterium}{
Information criterion used to select the number of components of the ME fit. One of \code{"AIC"} and \code{"BIC"} (default).
}
}
\details{
See Beirlant et al. (2015) and Verbelen et al. (2014) for details. The code follows the notation of the latter.
}
\value{
A \code{\link{SpliceFit}} object.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for fitting the mixed Erlang distribution.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}}, \code{\link[=pSplice]{Splice}},
\code{\link{Hill}}, \code{\link{trHill}}
}
\examples{
\dontrun{


# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.6)



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
