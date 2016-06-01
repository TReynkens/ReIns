\name{Hill.kopt}
\alias{Hill.kopt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select optimal threshold for Hill estimator
}
\description{
Select optimal threshold for the Hill estimator by minimising the Asymptotic Mean Squared Error (AMSE). 
}
\usage{
Hill.kopt(data, start = c(1, 1, 1), warnings = FALSE, logk = FALSE, 
          plot = FALSE, add = FALSE, main = "AMSE plot", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Vector of \eqn{n} observations. 
}
  \item{start}{
A vector of length 3 containing starting values for the first numerical optimisation (see \code{\link{Hill.2oQV}} for more details). Default is \code{c(1,1,1)}.
}
  \item{warnings}{
Logical indicating if possible warnings from the optimisation function are shown, default is \code{FALSE}.
}
  \item{logk}{
Logical indicating if the AMSE values are plotted as a function of the \eqn{log(k)} (\code{logk=TRUE}) or as a function of \eqn{k}. Default is \code{FALSE}.
}
  \item{plot}{
Logical indicating if the AMSE values should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the optimal value for \eqn{k} should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"AMSE plot"}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
See Section 4.2.1 of Albrecher et al. (2016) for more details.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{AMSE}{Vector of the AMSE values for each value of \eqn{k}.}
\item{kopt}{Optimal value of \eqn{k} corresponding to minimal \eqn{AMSE} value.}
\item{gammaopt}{Optimal value of the Hill estimator corresponding to minimal \eqn{AMSE} value.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Beirlant J., Goegebeur Y., Segers, J. and Teugels, J. (2004). \emph{Statistics of Extremes: Theory and Applications}, Wiley Series in Probability, Wiley, Chichester.

Beirlant J., Vynckier, P. and Teugels, J. (1996). "Tail Index Estimation, Pareto Quantile Plots, and Regression Diagnostics." \emph{Journal of the American Statistical Association}, 91, 1659--1667.
}
\author{
Tom Reynkens based on \code{S-Plus} code from Yuri Goegebeur.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Hill}}, \code{\link{Hill.2oQV}}
}
\examples{
data(norwegianfire)

# Plot Hill estimator as a function of k
Hill(norwegianfire$size[norwegianfire$year==76],plot=TRUE)

# Add optimal value of k
Hill.kopt(norwegianfire$size[norwegianfire$year==76],add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

