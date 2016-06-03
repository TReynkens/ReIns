\name{secura}
\alias{secura}
\docType{data}
\title{
Secura dataset
}
\description{
Automobile claims from 1988 to 2001, gathered from several European insurance companies, exceeding 1 200 000 Euro. Note that the data were, among others, corrected for inflation.
}
\usage{data("secura")}
\format{
  A data frame with 371 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{Year of claim occurence.}
    \item{\code{size}}{Size of automobile insurance claim (in EUR).}
  }
}

\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\examples{
data(secura)

# Exponential QQ-plot of Secura data
ExpQQ(secura$size)

# Pareto QQ-plot of Secura data
ParetoQQ(secura$size)

#Mean excess plot of Secura data (function of k)
MeanExcess(secura$size, k=TRUE)

# Mean excess plot of Secura data (function of order statistics)
MeanExcess(secura$size, k=FALSE)
}
\keyword{datasets}
