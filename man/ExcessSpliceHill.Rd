\name{ExcessSpliceHill}
\alias{ExcessSpliceHill}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using splicing of mixed Erlang and Pareto
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{u} using splicing of mixed Erlang and (truncated) Pareto distributions.
}
\usage{
ExcessSpliceHill(u, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
The retention level of the (re-)insurance or a vector of retention levels for the (re-)insurance.
}
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}.
}
}
\details{
This function is does also work for double splicing as described in \code{\link{SpliceFitHill}}.

See Beirlant et al. (2015) and Verbelen et al. (2014) for details.
}
\value{
An estimate for the premium is returned (for every value of \code{u}).
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
\code{\link{SpliceFitHill}}, \code{\link{ExcessHill}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.8)

# Excess-loss premium 
ExcessSpliceHill(2, splicefit)

}
}