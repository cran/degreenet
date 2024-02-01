%  File degreenet/man/reedmolloy.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{reedmolloy}
\alias{reedmolloy}
\title{Generate a undirected network with a given sequence of degrees}
\description{
Generate a undirected network where the degree of each actor is specified. The degree is the number of actors the actor is tied to.
This returns a \code{network} object and requires the \code{igraph} package.
}
\usage{
reedmolloy(deg, maxit=10, verbose=TRUE)
}
\arguments{
\item{deg}{vector of counts where element \eqn{i} is the degree of actor \eqn{i}. Its sum should be even.}
\item{maxit}{integer; maximum number of jitterings of the degree sequence to find a valid network.}
\item{verbose}{Print out details of the progress of the algorithm.}
}
\value{
 The network is returned as a \code{network} object.
}
\seealso{ayulemle, dyule}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for
details}
\examples{
# Now, simulate a Poisson Lognormal distribution over 100
# observations with mean = -1 and s.d. = 1.

set.seed(2)
s4 <- simpln(n=100, v=c(-1,1))
table(s4)
#
simr <- reedmolloy(s4)
simr
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
