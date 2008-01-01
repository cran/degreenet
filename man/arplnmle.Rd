\name{rplnmle}
\alias{rplnmle}
\alias{drpln}
\alias{llrplnall}
\alias{bootstraprpln}
\alias{ldrpln}
\alias{llrpln}
\alias{llrplnall}
\alias{rplnmle}
\title{Rounded Poisson Lognormal Modeling of Discrete Data}
\description{
Functions to Estimate the Rounded Poisson Lognormal Discrete Probability Distribution via maximum likelihood.
}
\usage{
rplnmle(x, cutoff = 1, cutabove = 1000, guess = c(-1,1),
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
\seealso{aplnmle}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
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

s4est <- rplnmle(s4)
s4est

}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
