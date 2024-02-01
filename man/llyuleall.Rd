%  File degreenet/man/llyuleall.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{llyuleall}
\alias{llyuleall}
\alias{lldpall}
\alias{llgeoall}
\alias{llgpall}
\alias{llgwall}
\alias{llgy0all}
\alias{llgyall}
\alias{llnb0all}
\alias{llnball}
\alias{llnbwall}
\alias{llnby0all}
\alias{llnbyall}
\alias{llnbzeroall}
\alias{llpeall}
\alias{llpoiall}
\alias{llsgeoall}
\title{Calculate the log-likelihood for Count Distributions}
\description{
Functions to Estimate the Log-likelihood for Discrete Probability Distributions.
}
\usage{
llyuleall(v, x, cutoff = 2, cutabove = 1000,  np=1)
}
\arguments{
\item{v}{A vector of parameters for the Yule (a 1-vector - the scaling exponent).}
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{np}{wnumber of parameters in the model. For the Yule this is 1.}
}
\value{
 the log-likelihood for the data \code{x} at parameter value \code{v}.
}
\seealso{ayulemle, llyule, dyule, llwarall}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
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

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model (i.e., rho=4, p=2/3)
#

s4warest <- awarmle(s4)
s4warest

#
# Compare the AICC and BIC for the two models
#

llyuleall(v=s4est$theta,x=s4)
llwarall(v=s4warest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
