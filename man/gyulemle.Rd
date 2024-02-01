%  File degreenet/man/gyulemle.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{gyulemle}
\alias{gyulemle}
\alias{gwarmle}
\title{Models for Count Distributions}
\description{
Functions to Estimate Parametric Discrete Probability Distributions via maximum likelihood
Based on categorical response
}
\usage{
gyulemle(x, cutoff = 1, cutabove = 1000, guess = 3.5, conc = FALSE, 
    method = "BFGS", hellinger = FALSE, hessian=TRUE)
}
\arguments{
\item{x}{A vector of categories for counts (one per observation). The values of \code{x}
and the categories are: \code{0=0, 1=1, 2=2, 3=3, 4=4, 5=5-10, 6=11-20, 7=21-100, 8=>100}}
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
 \item{result}{vector of parameter estimates - lower 95\% confidence value, upper 95\% confidence value, the PDF MLE, the asymptotic standard error,
and the number of data values >=cutoff and <=cutabove.}
 \item{theta}{The Yule MLE of the PDF exponent.}
 \item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
%\seealso{gyulemle, gwarmle, dyule, }
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

#
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
