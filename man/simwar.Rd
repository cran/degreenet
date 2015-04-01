%  File degreenet/man/simwar.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{simwar}
\alias{simwar}
\title{Simulate from a Waring Distribution}
\description{
Functions to generate random samples from a Waring Probability Distribution
}
\usage{
simwar(n=100, v=c(3.5, 0.1), maxdeg=10000)
}
\arguments{
\item{n}{number of samples to draw.}
\item{v}{Waring parameters: scaling exponent and probability of a new actor.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution).}
}
\value{
 vector of random draws or samples.
}
\seealso{awarmle, dwar}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Waring distribution over 100
# observations with Waring with exponent 3.5 and probability of a new
# actor 0.1.

set.seed(1)
s4 <- simwar(n=100, v=c(3.5, 0.1))
table(s4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
