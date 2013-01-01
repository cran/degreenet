%  File degreenet/man/acmpmle.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{acmpmle}
\alias{acmpmle}
\alias{dcmp}
\alias{llcmpall}
\alias{cmp.mutonatural}
\alias{cmp.naturaltomu}
\alias{dcmp.natural}
\alias{dcmp_mu}
\alias{ldcmp.natural}
\alias{llcmp}
\alias{llcmpall}
\title{Conway Maxwell Poisson Modeling of Discrete Data}
\description{
Functions to Estimate the Conway Maxwell Poisson Discrete Probability Distribution via maximum likelihood.
}
\usage{
acmpmle(x, cutoff = 1, cutabove = 1000, guess=c(7,3),
    method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE)
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
}
\value{
 \item{theta}{vector of MLE of the parameters.}
 \item{asycov}{asymptotic covariance matrix.}
 \item{asycor}{asymptotic correlation matrix.}
 \item{se}{vector of standard errors for the MLE.}
%\item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
\seealso{ayulemle, awarmle, simcmp}
\note{
See the working papers on \url{http://www.csss.washington.edu/Papers} for
details.

Based on the C code in the package \pkg{compoisson} written by
Jeffrey Dunn (2008).}
\examples{

# Simulate a Conway Maxwell Poisson distribution over 100
# observations with mean of 7 and variance of 3
# This leads to a mean of 1

set.seed(1)
s4 <- simcmp(n=100, v=c(7,3))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameters
#

acmpmle(s4)

}
\references{
{compoisson: Conway-Maxwell-Poisson Distribution},
{Jeffrey Dunn}, {2008}, {R package version 0.3}
}
\keyword{models}
