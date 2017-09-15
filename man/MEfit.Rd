\name{MEfit}
\alias{MEfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixed Erlang fit
}
\description{
Create an S3 object using a Mixed Erlang (ME) fit.
}
\usage{
MEfit(p, shape, theta, M, M_initial = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Vector of mixing weights, denoted by \eqn{\alpha} in Verbelen et al. (2015).
}
  \item{shape}{
Vector of shape parameters \eqn{r}.
}
  \item{theta}{
Scale parameter \eqn{\theta}.
}
  \item{M}{
Number of mixture components.
}
  \item{M_initial}{
Initial value provided for \code{M}. When \code{NULL} (default), not included in the object.
}
}
\details{
The rate parameter \eqn{\lambda} used in Albrecher et al. (2017) is equal to \eqn{1/\theta}.

See Reynkens et al. (2017) and Section 4.3 of Albrecher et al. (2017) for more details
}
\value{
An S3 object which contains the input arguments in a list.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2017). "Modelling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{Insurance: Mathematics and Economics}, 77, 65--77.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFit}}, \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}}, \code{\link{SpliceFitGPD}}
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

