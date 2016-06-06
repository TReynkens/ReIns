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
Output list from \code{\link{SpliceFitHill}}, \code{\link{SpliceFiticPareto}} or \code{\link{SpliceFitGPD}}.
}
}
\details{ 
The premium for the excess-loss insurance with retentation \eqn{R} and limit \eqn{L} is given by
\deqn{E(\min{(X-R)_+, L}) = I(R) - I(R+L)}
where \eqn{I(u)=\int_u^{\infty} (1-F(z)) dz} is the integrated tail function.  When \eqn{L=\infty}, the premium is equal to \eqn{I(R)}.

See Section 3.3 of Albrecher et al. (2016) for details.

}
\value{
An estimate for the premium is returned (for every value of \code{R}).
}
\references{
Albrecher, H., Beirlant, J. and Teugels, J. (2016). \emph{Reinsurance: Actuarial and Statistical Aspects}, Wiley, Chichester.

Verbelen, R., Gong, L., Antonio, K., Badescu, A. and Lin, S. (2015). "Fitting Mixtures of
Erlangs to Censored and Truncated Data Using the EM Algorithm." \emph{Astin Bulletin}, 45, 729--758
}
\author{
Tom Reynkens with \code{R} code from Roel Verbelen for the excess-loss of the mixed Erlang distribution.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SpliceFitHill}}, \code{\link{SpliceFiticPareto}},  \code{\link{SpliceFitGPD}}  
}
\examples{
\dontrun{

# Pareto random sample
X <- rpareto(1000, shape = 2)

# Splice ME and Pareto
splicefit <- SpliceFitHill(X, 0.8)

# Excess-loss premium 
ExcessSplice(R=2, L=Inf, splicefit=splicefit)

}
}

