\name{SpliceFit}
\alias{SpliceFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Splicing fit
}
\description{
Create an S3 object using \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitPareto}} or
\code{\link{SpliceFitGPD}}.
}
\usage{
SpliceFit(const, trunclower, t, type, MEfit, EVTfit)
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
Vector of types of the distributions:  \code{"ME"} and then for each fitted EVT distribution: \code{Pa} (Pareto), \code{tPa} (truncated Pareto), \code{"cPa"} (Pareto with right censoring), \code{"ciPa"} (Pareto with interval censoring), \code{"tciPa"} (Pareto with truncation and interval censoring) or \code{GPD} (GPD).
}
\item{MEfit}{
\code{\link{MEfit}} object with details on the mixed Erlang fit.
}
\item{EVTfit}{
\code{\link{EVTfit}} object with details on the EVT fit.
}
}
\details{
See Albrecher et al. (2015) for more details.
}
\value{
An S3 object containg the above input arguments and values for \eqn{\pi}, the splicing weights.
These splicing weights are equal to \deqn{\pi_1=const_1,   \pi_2=const_2-const_1, ...,\pi_{l+1}=1-const_l=1-(\pi_1+...+\pi_l)} when \eqn{l\ge 2} and
\deqn{\pi_1=const_1, \pi_2=1-const_1=1-\pi_1} when \eqn{l=1} where  \eqn{l} is the length of \code{const}.

A summary method is available.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2015). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2014). "Fitting Mixtures of Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{FEB Research Report KBI_1401A}.
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MEfit}}, \code{\link{EVTfit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFitcHill}}, \code{\link{SpliceFitGPD}}
}
\examples{
# Create MEfit object
mefit <- MEfit(p=c(0.65,0.35), shape=c(39,58), theta=16.19, M=2)

# Create EVTfit object
evtfit <- EVTfit(gamma=c(0.76,0.64))

# Create SpliceFit object
splicefit <- SpliceFit(const=c(0.5,0.996), trunclower=0, t=c(1020,39096), type=c("ME","Pa","Pa"),
                       MEfit=mefit, EVTfit=evtfit)

# Show summary
summary(splicefit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}
