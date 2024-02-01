%  File degreenet/man/bsdp.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{bsdp}
\alias{bsdp}
\alias{bootstrapdp}
\title{Calculate Bootstrap Estimates and Confidence Intervals for the Discrete
Pareto Distribution}
\description{
Uses the parametric bootstrap to estimate the bias and confidence
interval of the MLE of the Discrete Pareto Distribution.
}
\usage{
bsdp(x, cutoff=1, m=200, np=1, alpha=0.95)
bootstrapdp(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="adpmle")
}
\arguments{
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{m}{Number of bootstrap samples to draw.}
\item{np}{Number of parameters in the model (1 by default).}
\item{alpha}{Type I error for the confidence interval.}
\item{hellinger}{Minimize Hellinger distance of the parametric model from the
      data instead of maximizing the likelihood.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{guess}{Initial estimate at the MLE.}
\item{mle.meth}{Method to use to compute the MLE.}
}
\value{
 \item{dist}{matrix of sample CDFs, one per row.}
 \item{obsmle}{The Discrete Pareto MLE of the PDF exponent.}
 \item{bsmles}{Vector of bootstrap MLE.}
 \item{quantiles}{Quantiles of the bootstrap MLEs.}
 \item{pvalue}{p-value of the Anderson-Darling statistics relative to the bootstrap MLEs.}
 \item{obsmands}{Observed Anderson-Darling Statistic.}
 \item{meanmles}{Mean of the bootstrap MLEs.}
 \item{guess}{Initial estimate at the MLE.}
 \item{mle.meth}{Method to use to compute the MLE.}
}
\seealso{anbmle, simdp, lldp}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{
\dontrun{
# Now, simulate a Discrete Pareto distribution over 100
# observations with expected count 1 and probability of another
# of 0.2

set.seed(1)
s4 <- simdp(n=100, v=3.31)
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameter.
#

s4est <- adpmle(s4)
s4est

#
# Use the bootstrap to compute a confidence interval rather than using the 
# asymptotic confidence interval for the parameter.
#

bsdp(s4, m=20)
}
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
