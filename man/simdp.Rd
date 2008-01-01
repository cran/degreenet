\name{simdp}
\alias{simdp}
\title{Simulate from a Discrete Pareto Distribution}
\description{
Functions to generate random samples from a Discrete Pareto Probability Distribution
}
\usage{
simdp(n=100, v=3.5, maxdeg=10000)
}
\arguments{
\item{n}{number of samples to draw.}
\item{v}{Discrete Pareto parameters: PDF exponent.}
\item{maxdeg}{Maximum degree to sample (using truncation of the distribution).}
}
\value{
 vector of random draws or samples.
}
\seealso{adpmle, ddp}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Discrete Pareto distribution over 100
# observations with lognormal mean -1 and lognormal standard deviation 1.

set.seed(1)
s4 <- simdp(n=100, v=3.5)
table(s4)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
