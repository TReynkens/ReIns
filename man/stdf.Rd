\name{stdf}
\alias{stdf}
\alias{stdf2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-parametric estimators of the STDF
}
\description{
Non-parametric estimators of the stable tail dependence function (STDF): \eqn{\hat{l}_k(x)} and \eqn{\tilde{l}_k(x)}.
}
\usage{
stdf(x, k, X, alpha = 0.5)

stdf2(x, k, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \eqn{d}-dimensional point to estimate the STDF in.
}
  \item{k}{
Value of the tail index \eqn{k}.
}
  \item{X}{
A data matrix of dimensions \eqn{n} by \eqn{d} with observations in the rows.
}
  \item{alpha}{
The parameter \eqn{\alpha} of the estimator \eqn{\hat{l}_k(x)} (\code{stdf}), default is 0.5. This argument is not used in \code{stdf2}.
}
}
\details{
The stable tail dependence function in \eqn{x} can be estimated by
\deqn{ \hat{l}_k(x) = 1/k \sum_{i=1}^n 1_{\{\exists j\in\{1,\ldots, d\}: \hat{F}_j(X_{i,j})>1-k/n x_j\}}}
with 
\deqn{\hat{F}_j(X_{i,j})=(R_{i,j}-\alpha)/n}
where \eqn{R_{i,j}} is the rank of \eqn{X_{i,j}} among the \eqn{n} observations in the \eqn{j}-th dimension:
\deqn{R_{i,j}=\sum_{m=1}^n 1_{\{X_{m,j}\le X_{i,j}\}}.}
This estimator is implemented in \code{stdf}.

The second estimator is given by
\deqn{ \tilde{l}_k(x) = 1/k \sum_{i=1}^n 1_{\{X_{i,1}\ge X^{(1)}_{n-[kx_1]+1,n} or \ldots or X_{i,d}\ge X^{(d)}_{n-[kx_d]+1,n}\}}}
where \eqn{X_{i,n}^{(j)}} is the \eqn{i}-th smallest observation in the \eqn{j}-th dimension.
This estimator is implemented in \code{stdf2}.

See Section 4.5 of Beirlant et al. (2016) for more details.
}
\value{
\code{stdf} returns the estimate \eqn{\hat{l}_k(x)} and \code{stdf2} returns the estimate \eqn{\tilde{l}_k(x)}.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Generate data matrix
X <- cbind(rpareto(100,2), rpareto(100,3))

# Tail index
k <- 20

# Point to evaluate the STDF in
x <- c(2,3)

# First estimate
stdf(x, k, X)

# Second estimate
stdf2(x, k, X)
}
