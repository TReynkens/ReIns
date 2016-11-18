\name{SpliceFit}
\alias{SpliceFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing fit
}
\description{
Create an S3 object using \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or
\code{\link{SpliceFitGPD}}.
}
\usage{
SpliceFit(const, trunclower, t, type, MEfit, EVTfit, loglik = NULL, IC = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{const}{
Vector of splicing constants or a single splicing constant.
}
\item{trunclower}{
Lower truncation point.
}
\item{t}{
Vector of splicing point or a single splicing point.
}
\item{type}{
Vector of types of the distributions:  \code{"ME"} and then for each fitted EVT distribution: \code{"Pa"} (Pareto), \code{"TPa"} (truncated Pareto) or \code{"GPD"} (GPD).
}
\item{MEfit}{
\code{\link{MEfit}} object with details on the mixed Erlang fit.
}
\item{EVTfit}{
\code{\link{EVTfit}} object with details on the EVT fit.
}
\item{loglik}{
Log-likelihood of the fitted model. When \code{NULL} (default), not included in the object. 
}
\item{IC}{
Information criteria of the fitted model. When \code{NULL} (default), not included in the object. 
This vector should have length 1, 2 or 3 when included.
}
}
\details{
See Reynkens et al. (2016) and Section 4.3 in Albrecher et al. (2017) for details.
}
\value{
An S3 object containing the above input arguments and values for \eqn{\pi}, the splicing weights.
These splicing weights are equal to \deqn{\pi_1=const_1,   \pi_2=const_2-const_1, ...,\pi_{l+1}=1-const_l=1-(\pi_1+...+\pi_l)} when \eqn{l\ge 2} and
\deqn{\pi_1=const_1, \pi_2=1-const_1=1-\pi_1} when \eqn{l=1} where  \eqn{l} is the length of \code{const}.

A summary method is available.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2016). "Modeling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MEfit}}, \code{\link{EVTfit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceFitGPD}}
}
\examples{
# Create MEfit object
mefit <- MEfit(p=c(0.65,0.35), shape=c(39,58), theta=16.19, M=2)

# Create EVTfit object
evtfit <- EVTfit(gamma=c(0.76,0.64), endpoint=c(39096, Inf))

# Create SpliceFit object
splicefit <- SpliceFit(const=c(0.5,0.996), trunclower=0, t=c(1020,39096), type=c("ME","TPa","Pa"),
                       MEfit=mefit, EVTfit=evtfit)

# Show summary
summary(splicefit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
