\name{rSplice}
\alias{rSplice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random numbers from splicing distribution
}
\description{
Generate random numbers from the fitted splicing distribution.
}
\usage{
rSplice(n, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations.
}
  \item{splicefit}{
Output list from \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}} or \code{\link{SpliceFitGPD}}.
}
}
\details{ 
See Section 3.3 of Albrecher et al. (2015) for details.
}
\value{
A random sample of length \eqn{n}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qSplice}}, \code{\link{pSplice}}, \code{\link{dSplice}}, \code{\link{ES}}, \code{\link{SpliceFitHill}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}}
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.6)

# Random sample from splicing distribution
x <- rSplice(1000, splicefit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
