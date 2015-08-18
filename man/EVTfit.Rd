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
Vector of endpoints. When \code{NULL} (default), a vector containing \code{Inf} for each value of \code{gamma} will be used.
}
  \item{sigma}{
Vector of scale estimates for the GPD. When \code{NULL} (default), not included in the object. 
}
}
\details{
See Albrecher et al. (2015) for more details.
}
\value{
An S3 object containg the above input arguments.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFit}}, \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}}
}
\examples{
# Create MEfit object
mefit <- MEfit(p=c(0.65,0.35), shape=c(39,58), theta=16.19, M=2)

# Create EVTfit object
evtfit <- EVTfit(gamma=c(0.76,0.64))

# Create SpliceFit object
splicefit <- SpliceFit(const=c(0.5,0.996), trunclower=0, t=c(1020,39096), type=c("ME","Pa","Pa"),
                       MEfit=mefit, EVTfit=evtfit)

# Show summary
summary(splicefit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }

