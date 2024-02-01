%  File degreenet/man/degreenet-internal.Rd
%  Part of the statnet package, http://statnet.org
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
% Copyright 2003 Mark S. Handcock, University of California-Los Angeles
% Copyright 2007 The statnet Development Team
\name{degreenet-internal}
\alias{degreenet-internal}
\title{Internal degreenet Objects}
%
% These need to be documented!
%
\alias{adpmle}
\alias{adpmlef}
\alias{ageomle}
\alias{agpmle}
\alias{agwmle}
\alias{agy0mle}
\alias{agymle}
\alias{amle}
\alias{anb0mle}
\alias{anbmle}
\alias{anbwmle}
\alias{anby0mle}
\alias{anbymle}
\alias{apemle}
\alias{apoimle}
\alias{asgeomle}
\alias{ayulemlef}
\alias{reporting}
\alias{bootstrapcor}
\alias{bootstrapmix}
\alias{bootstrapnbconc}
\alias{bootstrappowexp}
\alias{bootstrapdpconc}
\alias{bootstrapgdp}
\alias{bootstrapgdpconc}
\alias{bootstrapgnbyconc}
\alias{bootstrapgyuleconc}
\alias{bootstrapnbyconc}
\alias{bootstrapyuleconc}
\alias{bootstraprdp}
\alias{bootstraprgw}
\alias{bootstraprwar}
\alias{bootstrapryule}
\alias{bootstrapwarconc}
\alias{bootstrapzipf}
\alias{bootstrapgnbconc}
\alias{bootstrapgwarconc}
\alias{bootstrapgyule}
\alias{dgeodp}
\alias{gdpmle}
\alias{ggeodpmle}
\alias{ggeomle}
\alias{ggymle}
\alias{gnbmle}
\alias{gnbymle}
\alias{gpoimle}
\alias{ldgeodp}
\alias{llggeo}
\alias{llggeoall}
\alias{llggeodp}
\alias{llggy}
\alias{llggyall}
\alias{llgnb}
\alias{llgnby}
\alias{llgnbyall}
%
\alias{dpln1}
\alias{dpln.refined}
\alias{dtp}
\alias{ldpln1}
\alias{ldtp}
%
\alias{llgplnall}
\alias{bootstrapgplnconc}
\alias{gplnmle}
\alias{llgpln}
\alias{llgpln}
%
\alias{aghdimle}
\alias{bootstrapgwar}
\alias{bootstrappe}
\alias{bootstrappeconc}
\alias{ddpe}
\alias{dghdi}
\alias{dnb}
\alias{dpe}
\alias{ldghdi}
\alias{llghdi}
\alias{llghdiall}
%
\alias{rcmp.mu}
%
% These are unlikely to need to be documented!
%
\alias{is.psd}
\alias{gauss.hermite}
\alias{hermite}
\alias{rmultinomial}
\alias{rmultz2}
\alias{ddp}
\alias{dgdp}
\alias{dgwar}
\alias{dgyule}
\alias{dgyule0}
\alias{dgyuleb}
\alias{dnbwar}
\alias{dnbyule}
\alias{dnbyule0}
\alias{dnbyuleb}
\alias{lddp}
\alias{ldgdp}
\alias{ldgwar}
\alias{ldgyule}
\alias{ldgyule0}
\alias{ldnbwar}
\alias{ldnbyule}
\alias{ldnbyule0}
\alias{llrdp}
\alias{llrdpall}
\alias{llrgw}
\alias{llrgwall}
\alias{llrgwf}
\alias{llrgwfall}
\alias{llrgwp}
\alias{llrgwpall}
\alias{llrgy}
\alias{llrgyall}
\alias{llrnb}
\alias{llrnball}
\alias{llrnbw}
\alias{llrnbwall}
\alias{llrwar}
\alias{llrwarall}
\alias{llryule}
\alias{llryuleall}
\alias{mands}
\alias{nbmean}
\alias{plotcdf}
\alias{polylog}
\alias{rdpmle}
\alias{reporting}
\alias{rgwfmle}
\alias{rgwmle}
\alias{rgwpmle}
\alias{rgymle}
\alias{rnbmle}
\alias{rnbwmle}
\alias{rwarmle}
\alias{ryulemle}
\alias{ryulemlef}
\alias{zeta}

\alias{NVL}
\alias{statnetStartupMessage}

\alias{compress.data.frame}
\alias{gcmpmle}
\alias{gllcmp}
\alias{gllcmp.narrow}
\alias{gllcmpall}
\alias{opttest}
\alias{sort.data.frame}

\description{
  Internal degreenet functions.
}
\usage{
adpmle(x, cutoff = 1, cutabove = 1000, guess = 3.5, hessian=TRUE)
}
\arguments{
\item{x}{A vector of counts (one per observation).}
\item{cutoff}{Calculate estimates conditional on exceeding this value.}
\item{cutabove}{Calculate estimates conditional on not exceeding this value.}
\item{guess}{Initial estimate at the MLE.}
\item{hessian}{Calculate the hessian of the information matrix (for use with
calculating the standard errors.}
}
\value{
 \item{theta}{vector of MLE of the parameters.}
 \item{asycov}{asymptotic covariance matrix.}
 \item{asycor}{asymptotic correlation matrix.}
 \item{se}{vector of standard errors for the MLE.}
%\item{value}{The maximized value of the function.}
 \item{conc}{The value of the concentration index (if calculated).}
}
\seealso{ayulemle, adqemle, simdqe}
\note{See the papers on \url{https://handcock.github.io/?q=Holland} for details.}
\references{Jones, J. H. and Handcock, M. S. "An assessment
of preferential attachment as a mechanism for human sexual
network formation," Proceedings of the Royal Society, B, 2003,
270, 1123-1128.}
\keyword{internal}
