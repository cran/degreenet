%  File degreenet/man/llpln.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{llpln}
\alias{llpln}
\title{Calculate the Conditional log-likelihood for the Poisson Lognormal Distributions}
\description{
Compute the Conditional Log-likelihood for the Poisson Lognormal
Discrete Probability Distribution. The likelihood is calculated
conditionl on the count being at least the cutoff value and less than
or equal to the cutabove value.
}
\usage{
llpln(v, x, cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE,logn = TRUE)
}
\arguments{
\item{v}{A vector of parameters for the Yule (a 1-vector - the scaling exponent).}
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{xr}{range of count values to use to approximate the set of all realistic counts.}
\item{hellinger}{Calculate the Hellinger distance of the parametric model 
from the data instead of the log-likelihood?}
\item{logn}{Use logn parametrization, that is, mean and variance on the
observation scale.}
}
\value{
the log-likelihood for the data \code{x} at parameter value \code{v}
(or the Hellinder distance if \code{hellinger=TRUE}).
}
\seealso{aplnmle, llplnall, dpln}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{

# Simulate a Poisson Lognormal distribution over 100
# observations with lognormal mean -1 and logormal standard deviation 1.

set.seed(1)
s4 <- simpln(n=100, v=c(-1,1))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for rho
#

s4est <- aplnmle(s4)
s4est

#
# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model
#

s4warest <- awarmle(s4)
s4warest

#
# Compare the log-likelihoods for the two models
#

llpln(v=s4est$theta,x=s4)
llwar(v=s4warest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
