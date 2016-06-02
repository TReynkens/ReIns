\name{pEdge}
\alias{pEdge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Edgeworth approximation

}
\description{
Edgeworth approximation of the CDF using the first four moments.
}
\usage{
pEdge(x, moments = c(0, 1, 0, 3), raw = TRUE, lower.tail = TRUE, log.p = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to approximate the CDF in.
}
  \item{moments}{
The first four raw moments if \code{raw=TRUE}, by default the first four raw moments of the standard normal distribution are used.
When \code{raw=FALSE}, the mean \eqn{\mu=E(X)}, variance \eqn{\sigma^2=E((X-\mu)^2)}, skewness (third standardised moment, \eqn{E((X-\mu)^3)/\sigma^3}) and kurtosis (fourth standardised moment, \eqn{E((X-\mu)^4)/\sigma^4}).
}
  \item{raw}{
When \code{TRUE} (default), the first four raw moments are provided in \code{moments}. Otherwise, the mean, variance, skewness and kurtosis are provided in \code{moments}. 
}
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as log(p), default is \code{FALSE}.
 }
}
\details{ 
Denote the standard normal PDF and CDF respectively by \eqn{\phi} and \eqn{\Phi}.
Let \eqn{\mu} be the first moment, \eqn{\sigma^2=E((X-\mu)^2)} the variance, \eqn{\mu_3=E((X-\mu)^3)} the third central moment and \eqn{\kappa_4=E((X-\mu)^4)} the fourth central moment of the random variable \eqn{X}.
The corresponding cumulants are given by \eqn{\kappa_1=\mu}, \eqn{\kappa_2=\sigma^2}, \eqn{\kappa_3=\mu_3} and \eqn{\kappa_4=\mu_4-3\sigma^4}. 


Now consider the random variable \eqn{Z=(X-\mu)/\sigma}, which has cumulants
0, 1, \eqn{\gamma_3=\kappa_3/\sigma^3}  and \eqn{\gamma_4=\kappa_4/\sigma^4=\mu_4/\sigma^4-3}.


The Edgeworth approximation for the CDF of \eqn{X} (\eqn{F(x)}) is given by
\deqn{\hat{F}_{E}(x) = \Phi(z) + \phi(z) (-\gamma_3/6 h_2(z)- (3\gamma_4h_3(z)+\gamma_3^2h_5(z))/72)}
with \eqn{h_2(z)=z^2-1}, \eqn{h_3(z)=z^3-3z}, \eqn{h_5(z)=z^5-10z^3+15z} and \eqn{z=(x-\mu)/\sigma}.


See Section 6.2 of Albrecher et al. (2016) for more details.
}
\value{
Vector of estimates for the probabilities \eqn{F(x)=P(X\le x)}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Cheah, P.K., Fraser, D.A.S. and Reid, N. (1993). "Some Alternatives to Edgeworth." \emph{The Canadian Journal of Statistics}, 21(2), 131--138.
}
\author{
Tom Reynkens
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pGC}}, \code{\link{pEdge}}
}
\examples{
# Chi-squared sample
X <- rchisq(1000, 2)


x <- seq(0, 10, 0.01)

# Empirical moments
moments = c(mean(X), mean(X^2), mean(X^3), mean(X^4))

# Gram-Charlier approximation
p1 <- pGC(x, moments)

# Edgeworth approximation
p2 <- pEdge(x, moments)

# Normal approximation
p3 <- pClas(x, mean(X), var(X))

# True probabilities
p <- pchisq(x, 2)


# Plot true and estimated probabilities
plot(x, p, type="l", ylab="F(x)", ylim=c(0,1), col="red")
lines(x, p1, lty=2)
lines(x, p2, lty=3)
lines(x, p3, lty=4, col="blue")

legend("bottomright", c("True CDF", "GC approximation", 
                        "Edgeworth approximation", "Normal approximation"), 
       col=c("red", "black", "black", "blue"), lty=1:4, lwd=2)
}