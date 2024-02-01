%  File degreenet/man/ryule.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{ryule}
\alias{ryule}
\title{Generate a (non-random) network from a Yule Distribution}
\description{
Generate a network with a given number of actors having a degree distribution draw
from a Yule distribution.
The resultant network is not random - that is, is not a random draw from all such networks.
}
\usage{
ryule(n=20,rho=2.5, maxdeg=n-1, maxit=10, verbose=FALSE)
}
\arguments{
\item{n}{Number of actors in the network.}
\item{rho}{PDF exponent of the Yule distribution.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution). If this is greater then \code{n-1} then \code{n-1} is used.}
\item{maxit}{integer; maximum number of resamplings of the degree sequence to find a valid network.}
\item{verbose}{Print out details of the progress of the algorithm.}
}
\value{
 If the \code{network} package is available, the network is returned as
a \code{network} object. If not a sociomatrix is returned.
}
\seealso{ayulemle, dyule, reedmolloy}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{
# Now, simulate a Yule network of 30
# actors with rho=4.0
ryule(n=30, rho=4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
