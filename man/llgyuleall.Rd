%  File degreenet/man/llgyuleall.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{llgyuleall}
\alias{llgyuleall}
\alias{llgwarall}
\alias{llgdpall}
\alias{llgnball}
\alias{llgnbwall}
\alias{llgpoiall}
\title{Calculate the log-likelihood for Count Distributions}
\description{
Functions to Estimate the Log-likelihood for Discrete Probability Distributions Based on Categorical Response.
}
\usage{
llgyuleall(v, x, cutoff = 2, cutabove = 1000,  np=1)
}
\arguments{
\item{v}{A vector of parameters for the Yule (a 1-vector - the scaling exponent).}
\item{x}{A vector of categories for counts (one per observation). The values of \code{x}
and the categories are: \code{0=0, 1=1, 2=2, 3=3, 4=4, 5=5-10, 6=11-20, 7=21-100, 8=>100}}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{np}{wnumber of parameters in the model. For the Yule this is 1.}
}
\value{
 the log-likelihood for the data \code{x} at parameter value \code{v}.
}
\seealso{gyulemle, llgyule, dyule, llgwarall}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{
#
# Simulate a Yule distribution over 100
# observations with rho=4.0
#

set.seed(1)
s4 <- simyule(n=100, rho=4)
table(s4)

#
# Recode it as categorical
#

s4[s4 >  4 & s4 < 11] <- 5
s4[s4 > 100] <- 8
s4[s4 >  20] <- 7
s4[s4 >  10] <- 6


#
# Calculate the MLE and an asymptotic confidence
# interval for rho
#

s4est <- gyulemle(s4)
s4est

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model (i.e., rho=4, p=2/3)
#

s4warest <- gwarmle(s4)
s4warest

#
# Compare the AICC and BIC for the two models
#

llgyuleall(v=s4est$theta,x=s4)
llgwarall(v=s4warest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
