\name{ExcessSplicecHill}
\alias{ExcessSplicecHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using splicing of mixed Erlang and Pareto adapted for censoring
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{u} using splicing of mixed Erlang and (truncated) Pareto distributions, adapted for censoring.
}
\usage{
ExcessSplicecHill(splicefit, u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}.
}
  \item{u}{
The retention level of the (re-)insurance.
}
}
\details{
See Beirlant et al. (2015) and Verbelen et al. (2014) for details.
}
\value{
An estimate for the premium is returned.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the excess-loss of the mixed Erlang distribution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFitcHill}}, \code{\link{ExcessHill}}, \code{\link{ExcessSpliceHill}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

Y <- rpareto(1000, shape = 3)

Z <- pmin(X,Y)

censored <- X>Y

# Splice ME and Pareto
splicefit <- SpliceFitcHill(Z, censored=censored, const=0.8)

# Excess-loss premium 
ExcessSplicecHill(splicefit, 1)

}
}