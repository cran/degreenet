\name{simpln}
\alias{simpln}
\title{Simulate from a Poisson Lognormal Distribution}
\description{
Functions to generate random samples from a Poisson Lognormal Probability Distribution
}
\usage{
simpln(n=100, v=c(-1,1), maxdeg=10000)
}
\arguments{
\item{n}{number of samples to draw.}
\item{v}{Poisson Lognormal parameters: lognormal mean and lognormal s.d.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution).}
}
\value{
 vector of random draws or samples.
}
\seealso{aplnmle, dpln}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Poisson Lognormal distribution over 100
# observations with lognormal mean -1 and lognormal standard deviation 1.

set.seed(1)
s4 <- simpln(n=100, v=c(-1,1))
table(s4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
