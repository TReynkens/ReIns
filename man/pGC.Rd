\name{pGC}
\alias{pGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gram-Charlier approximation
}
\description{
Gram-Charlier approximation of the CDF using the first four raw moments.
}
\usage{
pGC(x, moments = c(0, 1, 0, 3), lower.tail = TRUE, log.p = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to approximate the CDF in.
}
  \item{moments}{
The first four raw moments, by default the first four raw moments of the standard normal distribution are used.
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



The Gram-Charlier approximation for the CDF of \eqn{Z} (\eqn{F(z)}) is given by
\deqn{\hat{F}_{GC}(x) = \Phi(z) + \phi(z) (-\gamma_3/6 h_2(z)- \gamma_4/24h_3(z)}
with \eqn{h_2(z)=z^2-1} and \eqn{h_3(z)=z^3-3z}.

See Albrecher et al. (2016) for details.
}
\value{
Vector of estimates for the probabilities \eqn{F(x)=P(X\le x)}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Cheah, P.K., Fraser, D.A.S. and Reid, N. (1993). "Some Alternatives to Edgeworth." \emph{The Canadian Journal of Statistics}, 21(2), 131--138.
}
\author{
Tom Reynkens
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pEdge}}, \code{\link{pClas}}
}
\examples{
# Chi-squared sample
X <- rchisq(1000, 2)

# Empirical moments
moments = c(mean(X),mean(X^2),mean(X^3),mean(X^4))

x <- seq(0,10,0.01)

# Gram-Charlier approximation
p1 <- pGC(x, moments)

# Edgeworth approximation
p2 <- pEdge(x, moments)

# True probabilities
p <- pchisq(x, 2)


# Plot true and estimate probabilities
plot(x, p, type="l", ylab="F(x)", ylim=c(0,1))
lines(x, p1, lty=2)
lines(x, p2, lty=3)

legend("bottomright", c("True CDF", "GC approximation", 
                        "Edgeworth approximation"), lty=1:3)

}

