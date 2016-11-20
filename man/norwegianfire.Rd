\name{norwegianfire}
\alias{norwegianfire}
\docType{data}
\title{
Norwegian fire insurance data
}
\description{
Fire insurance claims for a Norwegian insurance company for the period 1972 to 1992 as studied in Beirlant et al. (1996). \cr A priority of 500 units was in force.

}
\usage{data("norwegianfire")}
\format{
  A data frame with 9181 observations on the following 2 variables.
  \describe{
    \item{\code{size}}{Size of fire insurance claim (in 1000 NOK).}
    \item{\code{year}}{Year of claim occurence (expressed as \eqn{yy} instead of \eqn{19yy}).}
  }
}

\source{
Beirlant, J., Teugels, J. L. and Vynckier, P. (1996). \emph{Practical Analysis of Extreme Values}, Leuven
University Press.
}
\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\examples{
data(norwegianfire)

# Exponential QQ-plot for Norwegian Fire Insurance data for claims in 1976.
ExpQQ(norwegianfire$size[norwegianfire$year==76])

# Pareto QQ-plot for Norwegian Fire Insurance data for claims in 1976.
ParetoQQ(norwegianfire$size[norwegianfire$year==76])
}
\keyword{datasets}
