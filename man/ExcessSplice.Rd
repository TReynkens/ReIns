\name{ExcessSplice}
\alias{ExcessSplice}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates for excess-loss premiums using splicing
}
\description{
Estimate premiums of excess-loss reinsurance with retention \eqn{R} and limit \eqn{L} using fitted spliced distribution.
}
\usage{
ExcessSplice(R, L, splicefit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
The retention level of the (re-)insurance or a vector of retention levels for the (re-)insurance.
}
  \item{L}{
The limit of the (re-)insurance or a vector of limits for the (re-)insurance.
}
  \item{splicefit}{
A \code{\link{SpliceFit}} object,  e.g. output from \code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
}
\details{ 
The premium for the excess-loss insurance with retention \eqn{R} and limit \eqn{L} is given by
\deqn{E(\min{(X-R)_+, L}) = \Pi(R) - \Pi(R+L)}
where \eqn{\Pi(u)=E((X-u)_+)=\int_u^{\infty} (1-F(z)) dz} is the premium of the excess-loss insurance with retention \eqn{u}.  When \eqn{L=\infty}, the premium is equal to \eqn{\Pi(R)}.

See Reynkens et al. (2016) and Section 4.5 of Albrecher et al. (2017) for more details.
}
\value{
An estimate for the premium is returned (for every value of \code{R}).
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2017). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Reynkens, T., Verbelen, R., Beirlant, J. and Antonio, K. (2016). "Modeling Censored Losses Using Splicing:
a Global Fit Strategy With Mixed Erlang and Extreme Value Distributions". \emph{arXiv:1608.01566}.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the estimates for the excess-loss premiums using the mixed Erlang distribution.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFitPareto}}, \code{\link{SpliceFiticPareto}},  \code{\link{SpliceFitGPD}}  
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitPareto(X, 0.8)

# Excess-loss premium 
ExcessSplice(R=2, L=Inf, splicefit=splicefit)
}
}

