%  File degreenet/man/aplnmle.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{aplnmle}
\alias{aplnmle}
\alias{dpln}
\alias{ldpln}
\alias{llplnall}
\title{Poisson Lognormal Modeling of Discrete Data}
\description{
Functions to Estimate the Poisson Lognormal Discrete Probability Distribution via maximum likelihood.
}
\usage{
aplnmle(x, cutoff = 1, cutabove = 1000, guess = c(0.6,1.2),
    method = "BFGS", conc = FALSE, hellinger = FALSE, hessian=TRUE,logn=TRUE)
}
\arguments{
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{guess}{Initial estimate at the MLE.}
\item{method}{Method of optimization. See "optim" for details.}
\item{conc}{Calculate the concentration index of the distribution?}
\item{hellinger}{Minimize Hellinger distance of the parametric model from the
data instead of maximizing the likelihood.}
\item{hessian}{Calculate the hessian of the information matrix (for use with
calculating the standard errors.}
\item{logn}{Use logn parametrization, that is, mean and variance on the observation scale.}
}
\value{
 \item{theta}{vector of MLE of the parameters.}
 \item{asycov}{asymptotic covariance matrix.}
 \item{asycor}{asymptotic correlation matrix.}
 \item{se}{vector of standard errors for the MLE.}
%\item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
\seealso{ayulemle, awarmle, simpln}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{

# Simulate a Poisson Lognormal distribution over 100
# observations with lognormal mean of -1 and lognormal variance of 1
# This leads to a mean of 1

set.seed(1)
s4 <- simpln(n=100, v=c(-1,1))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

s4est <- aplnmle(s4)
s4est

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Yule model
#

s4yuleest <- ayulemle(s4)
s4yuleest

# Calculate the MLE and an asymptotic confidence
# interval for rho under the Waring model
#

s4warest <- awarmle(s4)
s4warest

#
# Compare the AICC and BIC for the three models
#

llplnall(v=s4est$theta,x=s4)
llyuleall(v=s4yuleest$theta,x=s4)
llwarall(v=s4warest$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
