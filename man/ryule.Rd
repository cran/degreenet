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
ryule(n=20,rho=2.5,cutoff=1,cutabove=1000,
                  greedy=FALSE,
                  maxdeg=10000,maxout=TRUE,verbose=FALSE)
}
\arguments{
\item{n}{Number of actors in the network.}
\item{rho}{PDF exponent of the Yule distribution.}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{greedy}{Use a greedy algorithm that may perform better for larger networks,
and may perform worse.}
\item{maxdeg}{The maximum degree to simulate. The probabilities of the degrees are
chosen by truncation.}
\item{maxout}{Ensure the maximum degree of any actor is feasible - that is, at most
the number of available actors.}
\item{verbose}{Print out details of the progress of the algorithm.}
}
\value{
 If the \code{network} package is available, the network is returned as
a \code{network} object. If not a sociomatrix is returned.
}
\seealso{ayulemle, dyule, reedmolloy}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
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
