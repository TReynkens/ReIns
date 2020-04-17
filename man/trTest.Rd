\name{trTest}
\alias{trTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test for truncated Pareto-type tails
}
\description{
Test between non-truncated Pareto-type tails \emph{(light truncation)} and truncated Pareto-type tails \emph{(rough truncation)}.
}
\usage{
trTest(data, alpha = 0.05, plot = TRUE, main = "Test for truncation", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{alpha}{
 The used significance level, default is \code{0.05}.
}
 \item{plot}{
Logical indicating if the P-values should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Test for truncation"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics:plot.default]{plot}} for more details.
}
}
\details{
We want to test
\eqn{H_0: X} has non-truncated Pareto tails vs.
\eqn{H_1: X} has truncated Pareto tails.
Let
\deqn{E_{k,n}(\gamma) = 1/k \sum_{j=1}^k (X_{n-k,n}/X_{n-j+1,n})^{1/\gamma},}
with \eqn{X_{i,n}} the \eqn{i}-th order statistic.
The test statistic is then
\deqn{T_{k,n}=\sqrt{12k} (E_{k,n}(H_{k,n})-1/2) / (1-E_{k,n}(H_{k,n}))}
which is asymptotically standard normally distributed.
We reject \eqn{H_0} on level \eqn{\alpha} if 
\deqn{T_{k,n}<-z_{\alpha}}
where \eqn{z_{\alpha}} is the \eqn{100(1-\alpha)\%} quantile of a standard normal distribution.
The corresponding P-value is thus given by
\deqn{\Phi(T_{k,n})} with \eqn{\Phi} the CDF of a standard normal distribution.

See Beirlant et al. (2016) or Section 4.2.3 of Albrecher et al. (2017) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{testVal}{Corresponding test values.}
\item{critVal}{Critical value used for the test, i.e. \code{qnorm(1-alpha/2)}.}
\item{Pval}{Corresponding P-values.}
\item{Reject}{Logical vector indicating if the null hypothesis is rejected for a certain value of \code{k}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant, J., Fraga Alves, M.I. and Gomes, M.I. (2016). "Tail fitting for Truncated and Non-truncated Pareto-type Distributions." \emph{Extremes}, 19, 429--462.
}
\author{
Tom Reynkens.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trHill}}, \code{\link{trTestMLE}}
}
\examples{
# Sample from truncated Pareto distribution.
# truncated at 95\% quantile
shape <- 2
X <- rtpareto(n=1000, shape=shape, endpoint=qpareto(0.95, shape=shape))

# Test for truncation
trTest(X)


# Sample from truncated Pareto distribution.
# truncated at 99\% quantile
shape <- 2
X <- rtpareto(n=1000, shape=shape, endpoint=qpareto(0.99, shape=shape))

# Test for truncation
trTest(X)
}
\keyword{htest}