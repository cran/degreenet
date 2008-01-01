\name{simyule}
\alias{simyule}
\title{Simulate from a Yule Distribution}
\description{
Functions to generate random samples from a Yule Probability Distribution
}
\usage{
simyule(n=100, rho=4, maxdeg=10000)
}
\arguments{
\item{n}{number of samples to draw.}
\item{rho}{Yule PDF exponent.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution).}
}
\value{
 vector of random draws or samples.
}
\seealso{ayulemle, dyule}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Yule distribution over 100
# observations with rho=4.0

set.seed(1)
s4 <- simyule(n=100, rho=4)
table(s4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
