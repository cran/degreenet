#  File degreenet/R/cmp.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Calculate the CMP law MLE
#
acmpmle <-function(x,cutoff=1,cutabove=1000,guess=c(7,3),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- acmpmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 guess <- cmp.mutonatural(mu=guess[1],sig=guess[2])
 guess <- c(log(guess$lambda), log(guess$nu))
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llcmp,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llcmp,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  aaa$par <- cmp.naturaltomu(c(exp(aaa$par[1]),exp(aaa$par[2])))
  names(aaa$par) <- c("CMP mean","CMP s.d.")
  if(is.psd(-aaa$hessian)){
   asycov <- -solve(aaa$hessian)
   dimnames(asycov) <- list(names(aaa$par),names(aaa$par))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   names(asyse) <- names(aaa$par)
   ccc <- list(theta=aaa$par,asycov=asycov,se=asyse,asycor=asycor)
  }else{
   ccc <- list(theta=aaa$par)
  }
 }else{
  ccc <- list(theta=rep(NA,length=length(x)))
 }
#
 ccc$conc <- NA
 if(conc & exists("aaa")){
 v <- aaa$par
 probv <- v
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>0){
  c0 <- 1-sum(dcmp_mu(v=probv,x=0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dcmp_mu(v=probv,x=xr) / c0
 tx <- tabulate(x+1)/length(x)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nr <- tr[tr<cutoff]
 p0 <- sum(nc)
 p1 <- sum(nc * nr)
 p2 <- sum(nc * nr^2)
 conc<-(p1+sum(xr*pdf)*(1-p0))/(p2+sum(xr*xr*pdf)*(1-p0)-p1-sum(xr*pdf)*(1-p0))
 ccc$conc <- conc
 }
 ccc
}
#
# Complete data log-likelihoods
#
llcmpall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llcmp(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llcmp <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
 n <- length(x)
 if(cutoff>0){
  cprob <- 1 - sum(dcmp.natural(v=exp(probv),x=0:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- dcmp.natural(v=exp(probv),x=xr)
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldcmp.natural(v=exp(probv),x=xv))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
#
# Compute the CMP PMF
#
ldcmp.natural <- function(v,x,cutoff=0){
 log(dcmp.natural(v,x,cutoff=cutoff))
}
dcmp_mu <- function(v, x, cutoff=0, err=0.00001, log=FALSE){
  # Perform argument checking
  if (v[1] < 0 || v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  out <- cmp.mutonatural(v[1],v[2])
  y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(length(x)),
            ERr=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="degreenet")$val
  if (cutoff > 0) {
    c0 <- 1 - sum(dcmp_mu(v = v, x = (0:(cutoff - 1)), cutoff=0))
    y <- y/c0
  }
  return(y)
}
dcmp.natural <- function(v, x, cutoff=0, err=0.00001, log=FALSE){
  # Perform argument checking
  if (v[1] < 0 || v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(v[1]),
            nu=as.double(v[2]),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="degreenet")$val
  if (cutoff > 0) {
    c0 <- 1 - sum(dcmp.natural(v = v, x = (0:(cutoff - 1)), cutoff=0))
    y <- y/c0
  }
  return(y)
}
rcmp.mu <- function(n, mu, sig, err=0.00001, K=100){
  out <- cmp.mutonatural(mu,sig)
  # Perform argument checking
  if (out$lambda < 0 || out$nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="degreenet")
   return(out$x)
}
