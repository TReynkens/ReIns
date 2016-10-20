\name{pClas}
\alias{pClas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classical estimators for the CDF
}
\description{
Compute approximations of the CDF using the normal approximation,
normal-power approximation, shifted Gamma approximation or normal approximation to the shifted Gamma distribution.
}
\usage{
pClas(x, mean = 0, variance = 1, skewness = NULL, 
      method = c("normal", "normal-power", "shifted Gamma", "shifted Gamma normal"), 
      lower.tail = TRUE, log.p = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of points to approximate the CDF in.
}
  \item{mean}{
Mean of the distribution, default is 0.
}
  \item{variance}{
Variance of the distribution, default is 1.
}
  \item{skewness}{
Skewness coefficient of the distribution, this argument is not used for the normal approximation. Default is \code{NULL} meaning no skewness coefficient is provided.
}
  \item{method}{
Approximation method to use, one of \code{"normal"}, \code{"normal-power"}, \code{"shifted Gamma"} or \code{"shifted Gamma normal"}. Default is \code{"normal"}.
}
 \item{lower.tail}{
 Logical indicating if the probabilities are of the form \eqn{P(X\le x)} (\code{TRUE}) or \eqn{P(X>x)} (\code{FALSE}). Default is \code{TRUE.}
 }
  \item{log.p}{
  Logical indicating if the probabilities are given as \eqn{log(p)}, default is \code{FALSE}.
 }
}
\details{
\itemize{
\item The \emph{normal approximation} for the CDF of the r.v. \eqn{X} is defined as
\deqn{F_X(x) \approx \Phi((x-\mu)/\sigma)}
where \eqn{\mu} and \eqn{\sigma^2} are the mean and variance of \eqn{X}, respectively.


\item This approximation can be improved when the skewness parameter  
\deqn{\nu=E((X-\mu)^3)/\sigma^3}
is available. The \emph{normal-power approximation} of the CDF is then given by
\deqn{F_X(x) \approx \Phi(\sqrt{9/\nu^2 + 6z/\nu+1}-3/\nu)}
for \eqn{z=(x-\mu)/\sigma\ge 1} and \eqn{9/\nu^2 + 6z/\nu+1\ge 0}.

\item The \emph{shifted Gamma approximation} uses the approximation
\deqn{X \approx \Gamma(4/\nu^2, 2/(\nu*\sigma)) + \mu -2\sigma/\nu.}
Here, we need that \eqn{\nu>0}.

\item The \emph{normal approximation to the shifted Gamma distribution} approximates the CDF of \eqn{X} as
\deqn{F_X(x) \approx \Phi(\sqrt{16/\nu^2 + 8z/\nu}-\sqrt{16/\nu^2-1})}
for \eqn{z=(x-\mu)/\sigma\ge 1}. We need again that \eqn{\nu>0}.
}

See Section 6.2 of Albrecher et al. (2017) for more details.
}
\value{
Vector of estimates for the probabilities \eqn{F(x)=P(X\le x)}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.
}
\author{
Tom Reynkens
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pEdge}}, \code{\link{pGC}}
}
\examples{
# Chi-squared sample
X <- rchisq(1000, 2)

x <- seq(0, 10, 0.01)

# Classical approximations
p1 <- pClas(x, mean(X), var(X))
p2 <- pClas(x, mean(X), var(X), mean((X-mean(X))^3)/sd(X)^3, method="normal-power")
p3 <- pClas(x, mean(X), var(X), mean((X-mean(X))^3)/sd(X)^3, method="shifted Gamma")
p4 <- pClas(x, mean(X), var(X), mean((X-mean(X))^3)/sd(X)^3, method="shifted Gamma normal")

# True probabilities
p <- pchisq(x, 2)


# Plot true and estimated probabilities
plot(x, p, type="l", ylab="F(x)", ylim=c(0,1), col="red")
lines(x, p1, lty=2)
lines(x, p2, lty=3, col="green")
lines(x, p3, lty=4)
lines(x, p4, lty=5, col="blue")

legend("bottomright", c("True CDF", "normal approximation", "normal-power approximation",
                        "shifted Gamma approximation", "shifted Gamma normal approximation"), 
      lty=1:5, col=c("red", "black", "green", "black", "blue"), lwd=2)
}

