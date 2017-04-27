\name{EVTfit}
\alias{EVTfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EVT fit
}
\description{
Create an S3 object using an EVT (Extreme Value Theory) fit.
}
\usage{
EVTfit(gamma, endpoint = NULL, sigma = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
Vector of estimates for \eqn{\gamma}.
}
  \item{endpoint}{
Vector of endpoints (with the same length as \code{gamma}). When \code{NULL} (default), a vector containing \code{Inf} for each value of \code{gamma} will be used.
}
  \item{sigma}{
Vector of scale estimates for the GPD (with the same length as \code{gamma}). When \code{NULL} (default), not included in the object. 
}
}
\details{
See Reynkens et al. (2017) and Section 4.3 of Albrecher et al. (2017) for details.
}
\value{
An S3 object containing the above input arguments.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceFitGPD}}
}
\examples{
# Create MEfit object
mefit <- MEfit(p=c(0.65,0.35), shape=c(39,58), theta=16.19, M=2)

# Create EVTfit object
evtfit <- EVTfit(gamma=c(0.76,0.64), endpoint=c(39096, Inf))

# Create SpliceFit object
splicefit <- SpliceFit(const=c(0.5,0.996), trunclower=0, t=c(1020,39096), type=c("ME","TPa","Pa"),
                       MEfit=mefit, EVTfit=evtfit)

# Show summary
summary(splicefit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }

