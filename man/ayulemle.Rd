%  File degreenet/man/ayulemle.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{ayulemle}
\alias{ayulemle}
\alias{dyule}
\alias{ldyule}
\title{Yule Distribution Modeling of Discrete Data}
\description{
Functions to Estimate the Yule Discrete Probability Distribution via maximum likelihood.
}
\usage{
ayulemle(x, cutoff = 1, cutabove = 1000, guess = 3.5, conc = FALSE,
     method = "BFGS", hellinger = FALSE, hessian = TRUE, weights = rep(1, length(x)))
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
\item{weights}{sample weights on the observed counts.}
}
\value{
 \item{theta}{vector of MLE of the parameters.}
 \item{asycov}{asymptotic covariance matrix.}
 \item{asycor}{asymptotic correlation matrix.}
 \item{se}{vector of standard errors for the MLE.}
%\item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
\seealso{ayulemle, awarmle, simyule}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{

# Simulate a Yule distribution over 100
# observations with PDf exponent of 3.5

set.seed(1)
s4 <- simyule(n=100, rho=3.5)
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

s4est <- ayulemle(s4)
s4est

#
# Compute the AICC and BIC for the model
#

llyuleall(v=s4est$theta,x=s4)

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
