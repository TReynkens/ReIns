\name{soa}
\alias{soa}
\docType{data}
\title{
SOA Group Medical Insurance Large Claims Database
}
\description{
The SOA Group Medical Insurance Large Claims Database records, among others, all the claim amounts exceeding 25,000 USD in the year 1991.
}
\usage{data("soa")}
\format{
  A data frame with 75789 observations on the following variable:
  \describe{
    \item{\code{size}}{Claim size (in USD).}
  }
}

\source{
Grazier, K. L. and G'Sell Associates (1997). \emph{Group Medical Insurance Large Claims Database Collection and Analysis}. SOA Monograph M-HB97-1, Society of Actuaries, Schaumburg.

Society of Actuaries, \url{https://www.soa.org/resources/experience-studies/2000-2004/91-92-group-medical-claims/}.
}
\references{
Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\examples{
data(soa)

# Histogram of log-claim amount
hist(log(soa$size),breaks=seq(10,16,0.2),xlab="log(Claim size)")

# Exponential QQ-plot of claim amount
ExpQQ(soa$size)

# Mean excess plot of claim amount (function of k)
MeanExcess(soa$size, k=TRUE)

# Mean excess plot of claim amount (function of order statistics)
MeanExcess(soa$size, k=FALSE)
}
\keyword{datasets}
