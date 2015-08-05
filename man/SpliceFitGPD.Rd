\name{SpliceFitGPD}
\alias{SpliceFitGPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing of mixed Erlang and GPD using POT-MLE
}
\description{
Fit spliced distribution of a mixed Erlang and Generalised Pareto Distributions (GPD). The parameters of the GPD are determined using the POT-MLE approach. 
}
\usage{
SpliceFitGPD(X, const, M = 10, s = 1:10, trunclower = 0, ncores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data used for fitting the distribution.
}
  \item{const}{
Vector of length \eqn{l} containing the probabilities of the quantiles where the distributions will be spliced (splicing points). The ME distribution will be spliced with \eqn{l} GPD distributions.
}
  \item{M}{
Initial number of Erlang mixtures, default is 10. This number can change when determining
an optimal mixed Erlang fit using an AIC criterion.
}
  \item{s}{
Vector of spread factors for the EM algorithm, default is \code{1:10}. We loop over these factors when determining
an optimal mixed Erlang fit using an AIC criterion.
}
  \item{trunclower}{
Lower truncation point. Default is 0.
}
  \item{ncores}{
Number of cores to use when determining an optimal mixed Erlang fit using an AIC criterion.
When \code{NULL}, \code{max(nc-1,1)} cores are used where \code{nc} is the number of cores as determined by \code{\link{detectCores}}.
}
}
\details{
See Beirlant et al. (2015) and Verbelen et al. (2014) for details. The code follows the notation of the latter.
}
\value{
A list with components:
\item{MEfit}{List with details on the mixed Erlang fit.}
\item{EVTfit}{List with details on the EVT fit (GPD).}
\item{t}{Vector of splicing points.}
\item{trunclower}{Lower truncation point as given by the user.}
\item{const}{Vector of splicing constants, it can differ slightly from the input argument \code{const}.}
\item{type}{Vector containing the types of the distributions in the EVT fit: \code{"GPD"} for each GPD distribution that is used.}
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
\code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceCDF}}, \code{\link{SplicePDF}},
\code{\link{GPDfit}}
}
\examples{
\dontrun{

# GPD random sample
X <- rgpd(1000, gamma = 0.5, sigma = 2)

# Splice ME and GPD
splicefit <- SpliceFitGPD(X, 0.6)

x <- seq(0,20,0.01)
# Plot of splicing CDF
plot(x, SpliceCDF(x, splicefit), type="l", xlab="x", ylab="F(x)")
# Plot of splicing PDF
plot(x, SplicePDF(x, splicefit), type="l", xlab="x", ylab="f(x)")

# Fitted survival function and empirical survival function 
SpliceECDF(x, X, splicefit)

# Log-log plot with empirical survival function and fitted survival function
SpliceLL(x, X, splicefit)

# PP-plot of empirical survival function and fitted survival function
SplicePP(x, X, splicefit)

# PP-plot of empirical survival function and 
# fitted survival function with log-scales
SplicePP(x, X, splicefit,log=TRUE)

}
}
