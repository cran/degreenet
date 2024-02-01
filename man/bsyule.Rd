%  File degreenet/man/bsyule.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{bsyule}
\alias{bsyule}
\alias{bootstrapyule}
\title{Calculate Bootstrap Estimates and Confidence Intervals for the Yule
Distribution}
\description{
Uses the parametric bootstrap to estimate the bias and confidence
interval of the MLE of the Yule Distribution.
}
\usage{
bsyule(x, cutoff=1, m=200, np=1, alpha=0.95, v=NULL,
                   hellinger=FALSE, cutabove=1000)
bootstrapyule(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=3.31,hellinger=FALSE,
                          mle.meth="ayulemle")
}
\arguments{
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{m}{Number of bootstrap samples to draw.}
\item{np}{Number of parameters in the model (1 by default).}
\item{alpha}{Type I error for the confidence interval.}
\item{v}{Parameter value to use for the bootstrap distribution. By default it is the
         MLE of the data.}
\item{hellinger}{Minimize Hellinger distance of the parametric model from the
      data instead of maximizing the likelihood.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{guess}{Initial estimate at the MLE.}
\item{mle.meth}{Method to use to compute the MLE.}
}
\value{
 \item{dist}{matrix of sample CDFs, one per row.}
 \item{obsmle}{The Yule MLE of the PDF exponent.}
 \item{bsmles}{Vector of bootstrap MLE.}
 \item{quantiles}{Quantiles of the bootstrap MLEs.}
 \item{pvalue}{p-value of the Anderson-Darling statistics relative to the bootstrap MLEs.}
 \item{obsmands}{Observed Anderson-Darling Statistic.}
 \item{meanmles}{Mean of the bootstrap MLEs.}
}
\seealso{ayulemle, simyule, llyule}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{
# Now, simulate a Yule distribution over 100
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
# Use the bootstrap to compute a confidence interval rather than using the 
# asymptotic confidence interval for rho.
#

bsyule(s4, m=20)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
