\name{trParetoQQ}
\alias{trParetoQQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Truncated Pareto quantile plot
}
\description{
Extension of the Pareto QQ-plot as described in Beirlant et al. (2016).
}
\usage{
trParetoQQ(data, r = 1, DT, kstar = NULL, plot = TRUE, main = "tPa QQ-plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{r}{
Trimming parameter, default is \code{1} (no trimming).
}
\item{DT}{
Vector of \eqn{n-1} estimates for the odds ratio \eqn{D_T} obtained from \code{\link{trDT}}.
}
  \item{kstar}{
Value for \eqn{k} used to construct the plot. When \code{NULL} (default), a value will be chosen by maximising the
correlation between the empirical and theoretical quantiles (see Details).
}
  \item{plot}{
Logical indicating if the quantiles should be plotted in a Pareto QQ-plot, default is \code{TRUE}.
}
  \item{main}{
Title for the plot, default is \code{"tPa QQ-plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The Pareto QQ-plot for truncated data plots 
\deqn{ (-log(\hat{D}_{T,r,k^*,n}+j/(n+1)), log(X_{n-j+1,n}) )}
for \eqn{j=1,...,n-1}.

The value for \eqn{k^*} can be be given by the user or can be determined automatically.
In the latter case, we use the \eqn{k^*} that maximises the absolute value of the correlation between \eqn{log(X_{n-j+1,n})}
and \eqn{log(\hat{D}_{T,r,k^*,n}+j/(n+1))} for \eqn{j=1,...k} and \eqn{k^*>10}.

When taking \eqn{D_T=0}, one obtains the ordinary Pareto QQ-plot.

Note that the definition here differs slightly from the one in Beirlant et al. (2016).
We plot the empirical and theoretical quantiles the other way around and therefore have to add a minus (before the log).
}
\value{
A list with following components:
\item{pqq.the}{Vector of theoretical quantiles \eqn{-log(\hat{D}_{T,r,k^*,n}+j/(n+1))}, see Details.}
\item{pqq.emp}{Vector of the empirical quantiles from the log-transformed data.}
\item{kstar}{Optimal value for \eqn{k} or input argument \code{kstar}, see Details.}
\item{DTstar}{Estimate of \eqn{D_T} corresponding to \code{kstar}.}

}

\references{
Beirlant, J., Fraga Alves, M.I. and Gomes, M.I. (2016). "Tail fitting for Truncated
Pareto-type Distributions." \emph{Extremes}, accepted.
}
\author{
Tom Reynkens.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ParetoQQ}}, \code{\link{trDT}}
}
\examples{
# Endpoint of truncated Pareto distribution
endpoint <- qpareto(0.99,shape=2)

# Generate sample from truncated Pareto distribution
X <- rtpareto(1000,shape=2,endpoint=endpoint)

# Ordinary Pareto QQ-plot
ParetoQQ(X)

# Truncated Hill estimates
gamma <- trHill(X)$gamma

# Estimates for odds ratio
dt <- trDT(X, gamma=gamma)$DT

# Truncated Pareto QQ-plot
trParetoQQ(X, DT=dt)
}

