\name{reedmolloy}
\alias{reedmolloy}
\title{Generate a (non-random) network with a given degree sequence}
\description{
Generate a network with a given number of actors with a given degree. The resultant
network is not random - that is, is not a random draw from all such networks.
This requires the \code{network} package.
}
\usage{
reedmolloy(deg, greedy=FALSE, verbose=TRUE)
}
\arguments{
\item{deg}{vector of counts where element i is the number of actors with degree i.}
\item{greedy}{Use a greedy algorithm that may perform better for larger networks,
and may perform worse.}
\item{verbose}{Print out details of the progress of the algorithm.}
}
\value{
 If the \code{network} package is available, the network is returned as
a \code{network} object. If not a sociomatrix is returned.
}
\seealso{ayulemle, dyule}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Poisson Lognormal distribution over 100
# observations with mean = -1 and s.d. = 1.

set.seed(6)
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
