%  File degreenet/man/simnb.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{simnb}
\alias{simnb}
\title{Simulate from a Negative Binomial Distribution}
\description{
Functions to generate random samples from a Negative Binomial Probability Distribution
}
\usage{
simnb(n=100, v=c(5,0.2), maxdeg=10000)
}
\arguments{
\item{n}{number of samples to draw.}
\item{v}{Negative Binomial parameters: expected count and probability of another.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution).}
}
\value{
 vector of random draws or samples.
}
\seealso{anbmle, dnb}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Negative Binomial distribution over 100
# observations with lognormal mean -1 and lognormal standard deviation 1.

set.seed(1)
s4 <- simnb(n=100, v=c(5,0.2))
table(s4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
