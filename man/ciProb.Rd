\name{ciProb}
\alias{ciProb}
\alias{ciReturn}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimator of small exceedance probabilities and return periods using interval censored Hill
}
\description{
Computes estimates of a small exceedance probability \eqn{P(X>q)} or return period \eqn{1/P(X>q)} using the estimates for the EVI obtained from the Hill estimator adapted for interval censoring.
}
\usage{
ciProb(Z, I, censored, gamma1, q, plot = FALSE, add = FALSE, 
        main = "Estimates of small exceedance probability", ...)
        
ciReturn(Z, I, censored, gamma1, q, plot = FALSE, add = FALSE, 
        main = "Estimates of return period", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Z}{
Vector of length \eqn{n} with the left boundary of the intervals. In an 
insurance context, these are the paid amounts.
}
  \item{I}{
Vector of length \eqn{n} with the right boundary of the intervals. In an 
insurance context, these are the incurred values.
}
  \item{censored}{
 A logical vector of length \eqn{n} indicating if an observation is censored.
}
  \item{gamma1}{
Vector of \eqn{n-1} estimates for the EVI obtained from \code{\link{ciHill}}.
}
  \item{q}{
The used large quantile (we estimate \eqn{P(X>q)} or \eqn{1/P(X>q)} for \eqn{q} large).
}
 \item{plot}{
Logical indicating if the estimates should be plotted as a function of \eqn{k}, default is \code{FALSE}.
}
  \item{add}{
Logical indicating if the estimates should be added to an existing plot, default is \code{FALSE}.
}
  \item{main}{
Title for the plot, default is \code{"Estimates of small exceedance probability"} for \code{ciProb}
and \code{"Estimates of return period"} for \code{ciReturn}.
}
  \item{\dots}{
Additional arguments for the \code{plot} function, see \code{\link[graphics]{plot}} for more details.
}
}
\details{
The probability is estimated as
\deqn{ (1-tb) * (q/Z_{n-k,n})^{-1/H_{k,n}^c}}
with \eqn{Z_{i,n}} the \eqn{i}-th order statistic of the data, \eqn{H_{k,n}^c}
the Hill estimator adapted for right censoring and \eqn{tb} the Turnbull estimator for the CDF evaluated in \eqn{Z_{n-k,n}}.
}
\value{
A list with following components:
\item{k}{Vector of the values of the tail parameter \eqn{k}.}
\item{P}{Vector of the corresponding probability estimates, only returned for \code{ciProb}.}
\item{R}{Vector of the corresponding estimates for the return period, only returned for \code{ciReturn}.}
\item{q}{The used large quantile.}
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.
}
\author{
Tom Reynkens.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ciHill}}, \code{\link{cProb}}, \code{\link{Prob}}, \code{\link{Turnbull}}
}

\examples{
data(secura)

# Censoring indicator
censored <- sample(0:1,length(secura$size),replace=TRUE)

# Incurred values
I <- secura$size+rchisq(length(secura$size),2)

# Hill estimator adapted for interval censoring
cihill <- ciHill(secura$size,I,censored=censored,plot=TRUE)

# Small exceedance probability
q <- 10^7
ciProb(secura$size,I,censored=censored,gamma1=cihill$gamma1,q=q,plot=TRUE)     

# Return period
ciReturn(secura$size,I,censored=censored,gamma1=cihill$gamma1,q=q,plot=TRUE)    
}
