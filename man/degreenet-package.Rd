%  File degreenet/man/degreenet-package.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{degreenet-package}
\alias{degreenet-package}
\alias{degreenet}
\docType{package}
\title{
Models for Skewed Count Distributions Relevant to Networks
}
\description{
\pkg{degreenet} is a collection of functions to fit, diagnose, and simulate from distributions for skewed count data. The coverage of distributions is very selective, focusing on those that have been proposed to model the degree distribution on networks. For the rationale for this choice, see the papers in the references section below.
For a list of functions type: help(package='degreenet')

For a complete list of the functions, use \code{library(help="degreenet")} 
or read the rest of the manual. For a simple demonstration,
use \code{demo(packages="degreenet")}.

The \pkg{degreenet} package is part of the \pkg{statnet} suite of packages. The suite was developed to facilitate the statistical analysis of network data.

When publishing results obtained using this package alone see the 
citation in \code{citation(package="degreenet")}. The citation for the original
paper to use this package is Handcock and Jones (2003) and it should be cited
for the theoretical development.

If you use other packages in the \pkg{statnet} suite, please cite it
as:

Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau,
and Martina Morris. 2003 
\emph{statnet: Software tools for the Statistical Modeling of Network Data}
\cr
\url{https://statnet.org}.
For complete citation information, use\cr
\code{citation(package="statnet")}.

All programs derived from this or other \pkg{statnet} packages
must cite them appropriately.
}
\details{
See the Handcock and Jones (2003) reference (and the papers it cites and is
cited by) for more information on the methodology.

Recent advances in the statistical modeling of random networks have had an
impact on the empirical study of social networks. Statistical exponential
family models (Strauss and Ikeda 1990) are a generalization of the Markov
random network models introduced by Frank and Strauss (1986).
These models
recognize the complex dependencies within relational data structures.  To
date, the use of stochastic network models for networks has been limited by
three interrelated factors: the complexity of realistic models, the lack of
simulation tools for inference and validation, and a poor understanding of the
inferential properties of nontrivial models.

This package relies on the \code{\link[network]{network}} package which allows networks to be
represented in \R. The \pkg{statnet} suite of packages allows maximum likelihood estimates of
exponential random network models to be calculated using Markov Chain Monte
Carlo, as well as a broad range of statistical analysis of networks, such as
tools for plotting networks, simulating
networks and assessing model goodness-of-fit.

For detailed information on how to download and install the software,
go to the \pkg{statnet} website:
\url{https://statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Mark S. Handcock \email{handcock@stat.ucla.edu}

Maintainer: Mark S. Handcock \email{handcock@stat.ucla.edu}
}
\references{
Frank, O., and Strauss, D.(1986). \emph{Markov graphs.} Journal of the American
Statistical Association, 81, 832-842. 

Jones, J. H. and Handcock, M. S. (2003).
\emph{An assessment of preferential attachment as a mechanism for human sexual network formation}, Proceedings of the Royal Society, B, 2003, 270, 1123-1128.

Handcock, M. S., Hunter, D. R., Butts, C. T.,  Goodreau,
S. M., and Morris, M. (2003), 
\emph{statnet: Software tools for the Statistical Modeling of Network Data.},\cr
URL \url{https://statnet.org}

Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
networks. Journal of the American Statistical Association, 85, 204-212. 
}
\keyword{ package }
\keyword{ models }
