%  File degreenet/man/awarmle.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{awarmle}
\alias{awarmle}
\alias{dwar}
\alias{ldwar}
\alias{llwarall}
\title{Waring Modeling of Discrete Data}
\description{
Functions to Estimate the Waring Discrete Probability Distribution via maximum likelihood.
}
\usage{
awarmle(x, cutoff = 1, cutabove = 1000, guess = c(3.5,0.1),
    method = "BFGS", conc = FALSE, hellinger = FALSE, hessian=TRUE)
}
\arguments{
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{guess}{Initial estimate at the MLE.}
\item{conc}{Calculate the concentration index of the distribution?}
\item{method}{Method of optimization. See "optim" for details.}
\item{hellinger}{Minimize Hellinger distance of the parametric model from the
data instead of maximizing the likelihood.}
\item{hessian}{Calculate the hessian of the information matrix (for use with
calculating the standard errors.}
}
\value{
 \item{theta}{vector of MLE of the parameters.}
 \item{asycov}{asymptotic covariance matrix.}
 \item{asycor}{asymptotic correlation matrix.}
 \item{se}{vector of standard errors for the MLE.}
%\item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
\seealso{ayulemle, awarmle, simwar}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{

# Simulate a Waring distribution over 100
# observations with a PDf exponent of 3.5 and a 
# probability of including a new actor of 0.1

set.seed(1)
s4 <- simwar(n=100, v=c(3.5,0.1))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

s4est <- awarmle(s4)
s4est

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Yule model
#

s4yuleest <- ayulemle(s4)
s4yuleest

#
# Compare the AICC and BIC for the two models
#

llwarall(v=s4est$theta,x=s4)
llyuleall(v=s4yuleest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
