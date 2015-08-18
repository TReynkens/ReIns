\name{SpliceFitcHill}
\alias{SpliceFitcHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing of mixed Erlang and Pareto for censored data
}
\description{
Fit spliced distribution of a mixed Erlang and a Pareto distribution, both adapted for right or interval censoring. The shape parameter of the Pareto distribution is determined using the Hill estimator adpated for censoring. 
}
\usage{
SpliceFitcHill(Z, I = Z, censored, const, M = 3, s = 1:10, trunclower = 0,
              EVTtruncation = FALSE, ncores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{const}{
Probability of quantile where the distributions will be spliced (splicing point).
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
Logical indicating if a truncated (and censored) Pareto distribution should be fitted as second splicing distribution, default is \code{FALSE}. This is ignored for right-censored data!
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
\item{const}{The splicing constant, it can differ slightly from the input argument \code{const}.}
\item{trunclower}{Lower truncation point as given by the user.}
\item{t}{Splicing point.}
\item{type}{Type of the distributions, here \code{"ME"} and then \code{"cPa"} (right censoring), \code{"ciPa"} (interval censoring) or \code{"tciPa"} (truncation and interval censoring).}
\item{MEfit}{List with details on the mixed Erlang fit.}
\item{EVTfit}{List with details on the EVT fit (Pareto distribution).}
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
\code{\link{SpliceFitHill}}, \code{\link{SpliceFitGPD}}, \code{\link[=pSplice]{Splice}},
\code{\link{cHill}}, \code{\link{ciHill}}, \code{\link{trciHill}}
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

