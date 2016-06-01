\name{MeanExcess}
\alias{MeanExcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean excess function
}
\description{
Computes the mean excess values for a vector of observations. These mean excess values can then be plotted as a function of the data or as a function of the tail parameter \eqn{k}.
}
\usage{
MeanExcess(data, plot = TRUE, k = FALSE, main = "Mean excess plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Vector of \eqn{n} observations. 
}
  \item{plot}{
Logical indicating if the mean excess values should be plotted in a mean excess plot, default is \code{TRUE}.
}
  \item{k}{
Logical indicating if the mean excess scores are plotted as a function of the tail parameter \eqn{k} (\code{k=TRUE}) or as a function of the data (\code{k=FALSE}). Default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Mean excess plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The mean excess plot is \deqn{(k,e_{k,n})} or \deqn{(X_{n-k,n}, e_{k,n})}
with \deqn{e_{k,n}=1/k\sum_{j=1}^k X_{n-j+1,n}-X_{n-k,n}.}

See Section 4.1 of Albrecher et al. (2016) for more details. 
Note that the mean excess plot is the derivative plot of the Exponential QQ-plot.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \code{k}.}
\item{X}{Vector of the order statistics \code{data[n-k]} corresponding to the tail parameters in \code{k}.}
\item{e}{Vector of the mean excess values corresponding to the tail parameters in \code{k}.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}
\seealso{
\code{\link{ExpQQ}}, \code{\link{LognormalQQ_der}}, \code{\link{ParetoQQ_der}}, \code{\link{WeibullQQ_der}} 
}
\examples{
data(norwegianfire)

# Mean excess plots for Norwegian Fire Insurance data for claims in 1976.

# Mean excess values as a function of k
MeanExcess(norwegianfire$size[norwegianfire$year==76], k=TRUE)

# Mean excess values as a function of the data
MeanExcess(norwegianfire$size[norwegianfire$year==76], k=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

