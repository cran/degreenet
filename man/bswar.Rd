\name{bswar}
\alias{bswar}
\alias{bootstrapwar}
\title{Calculate Bootstrap Estimates and Confidence Intervals for the Waring Distribution}
\description{
Uses the parametric bootstrap to estimate the bias and confidence
interval of the MLE of the Waring Distribution.
}
\usage{
bswar(x, cutoff=1, m=200, np=2, alpha=0.95, v=NULL,
                   hellinger=FALSE)
bootstrapwar(x,cutoff=1,cutabove=1000,
             m=200,alpha=0.95,guess=c(3.31, 0.1),file="none",
             conc = FALSE)
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
\item{guess}{Guess at the parameter value.}
\item{file}{Name of the file to store the results. By default do not save the
            results.}
\item{conc}{Calculate the concentration index of the distribution?}
}
\value{
 \item{dist}{matrix of sample CDFs, one per row.}
 \item{obsmle}{The Waring MLE of the PDF exponent.}
 \item{bsmles}{Vector of bootstrap MLE.}
 \item{quantiles}{Quantiles of the bootstrap MLEs.}
 \item{pvalue}{p-value of the Anderson-Darling statistics relative to the bootstrap MLEs.}
 \item{obsmands}{Observed Anderson-Darling Statistic.}
 \item{meanmles}{Mean of the bootstrap MLEs.}
 \item{guess}{Initial estimate at the MLE.}
 \item{mle.meth}{Method to use to compute the MLE.}
}
\seealso{anbmle, simwar, llwar}
\note{See the working papers on \url{http://www.csss.washington.edu/Papers} for
details}
\examples{
# Now, simulate a Waring distribution over 100
# observations with expected count 1 and probability of another
# of 0.2

set.seed(1)
s4 <- simwar(n=100, v=c(5,0.2))
table(s4)

#
# Calculate the MLE and an asymptotic confidence
# interval for the parameter.
#

s4est <- awarmle(s4)
s4est

#
# Use the bootstrap to compute a confidence interval rather than using the 
# asymptotic confidence interval for the parameter.
#

bswar(s4, m=20)
}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{models}
