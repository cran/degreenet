%  File degreenet/man/llyule.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{llyule}
\alias{llyule}
\alias{llwar}
\alias{lldp}
\alias{lldp.good}
\alias{llgeo}
\alias{llgp}
\alias{llgp.good}
\alias{llgw}
\alias{llgy}
\alias{llgy0}
\alias{llnb}
\alias{llnb0}
\alias{llnbw}
\alias{llnby}
\alias{llnby0}
\alias{llnbzero}
\alias{llpe}
\alias{llpoi}
\alias{llsgeo}
\title{Calculate the Conditional log-likelihood for Count Distributions}
\description{
Functions to Estimate the Conditional Log-likelihood for Discrete Probability Distributions. The likelihood is calcualted condition on the count being at least the cutoff value and less than or equal to the cutabove value.
}
\usage{
llyule(v, x, cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE, weights = rep(1, length(x)))
}
\arguments{
\item{v}{A vector of parameters for the Yule (a 1-vector - the scaling exponent).}
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{xr}{range of count values to use to approximate the set of all realistic counts.}
\item{hellinger}{Calculate the Hellinger distance of the parametric model 
from the data instead of the log-likelihood?}
\item{weights}{sample weights on the observed counts.}
}
\value{
the log-likelihood for the data \code{x} at parameter value \code{v}
(or the Hellinder distance if \code{hellinger=TRUE}).
}
\seealso{ayulemle, llyuleall, dyule}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{

# Simulate a Yule distribution over 100
# observations with rho=4.0

set.seed(1)
s4 <- simyule(n=100, rho=4)
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for rho
#

s4est <- ayulemle(s4)
s4est

#
# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model (i.e., rho=4, p=2/3)
#

s4warest <- awarmle(s4)
s4warest

#
# Compare the log-likelihoods for the two models
#

llyule(v=s4est$theta,x=s4)
llwar(v=s4warest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
