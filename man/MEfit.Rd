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
Vector of mixing weights.
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
The rate paremeter \eqn{\lambda} used in Albrecher et al. (2016) is equal to \eqn{1/\theta}.

The mixing weights are denoted by \eqn{\alpha} (no truncation) or \eqn{\beta} (truncation) in Verbelen et al. (2015).
}
\value{
An S3 object which contains the input arguments in a list.
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{"Reinsurance: Actuarial and Statistical Aspects"}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFit}}, \code{\link{SpliceFitHill}}, \code{\link{SpliceFitciPareto}}, \code{\link{SpliceFitGPD}}
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

