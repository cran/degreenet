#
# Bootstrap CI for Poisson Lognormal
#
bootstrappln <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(-1,1),
                          file="none"){
mle <- aplnmle(x=x,cutoff=cutoff,guess=guess)$theta
bmles <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 bmles[i] <- aplnmle(x=xsamp,cutoff=cutoff,guess=guess)$theta[1]
}
#
if(!missing(file)){
 save(mle,bmles,file=file)
}
c(quantile(bmles,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),mle)
}
#
# Bootstrap CI for Poisson Lognormal Conc
#
bootstrapplnconc <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(-1,1),
                          file="none"){
mle <- aplnmle(x=x,cutoff=cutoff,guess=guess,conc=TRUE)
cmle <- mle$conc
mle <- mle$theta
bmles <- rep(0,length=m)
bconc <- rep(0,length=m)
for(i in seq(along=bmles)){
 xsamp <- sample(x,size=length(x),replace=TRUE)
 tbs <- aplnmle(x=xsamp,cutoff=cutoff,guess=mle,conc=TRUE)
 if(is.na(tbs$theta[1])){bmles[i] <- bmles[i-1]; bconc[i] <- bconc[i-1]}
 else{
  bmles[i] <- tbs$theta[1]
  bconc[i] <- tbs$conc
 }
}
#
if(!missing(file)){
 save(cmle,mle,bconc,bmles,file=file)
}
mle <- c(quantile(bconc,c(0.5*(1-alpha),(1-alpha),0.5,alpha,0.5+0.5*alpha),na.rm=TRUE),cmle,
  mean(bmles < 3.0,na.rm=TRUE))
names(mle)[7] <- "prob < 3"
mle
}
bspln <- function(x, cutoff=1, m=200, np=2, alpha=0.95, v=NULL,
                   hellinger=FALSE){
 if(missing(v)){v <- aplnmle(x=x,cutoff=cutoff,hellinger=hellinger)$theta}
 obsmands <- mands(x,type="pln",v=v,cutoff=cutoff,hellinger=hellinger)
 bsmles <- matrix(0,nrow=m,ncol=np+1)
 for(i in 1:m){
#
# Sample from a Poisson Lognormal CDF at the MLE
#
  xsamp <- sample(x=obsmands$cdf[,"k"], 
                  size=length(x),
                  replace=TRUE,
		  prob=diff(c(0,obsmands$cdf[,"cdf"]))
		 )
  wmle <- aplnmle(x=xsamp,cutoff=cutoff,hellinger=hellinger)$theta
  bsmles[i,-np-1] <- wmle
  bsmles[i, np+1] <- mands(xsamp,type="pln",v=wmle,
                      cutoff=cutoff,hellinger=hellinger)$ads
# if(done){print(paste("done",i))}
 }
 if(hellinger){
  dimnames(bsmles) <- list(1:m,c("LN mean","LN s.d.","Pen. Hell. Dist."))
 }else{
  dimnames(bsmles) <- list(1:m,c("LN mean","LN s.d.","MANDS"))
 }
 #
 list(dist=obsmands$cdf,
  obsmle=v,
  bsmles=bsmles,
  quantiles=quantile(bsmles[,np+1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),
  pvalue=sum(obsmands$ads < bsmles[,np+1],na.rm=TRUE)/(sum(bsmles[,np+1]>0,na.rm=TRUE)+1),
  obsmands=obsmands$ads,
  meanmles=apply(bsmles,2,mean)
 )
}
#
# Calculate the Poisson Lognormal law MLE
#
aplnmle <-function(x,cutoff=1,cutabove=1000,guess=c(-1,1), approxlim=10,
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- aplnmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llpln,
#  lower=1.1,upper=30,
#  method="L-BFGS-B",
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger,
   approxlim=approxlim)
  aaanm <- optim(par=guess,fn=llpln,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger,
   approxlim=approxlim)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Poisson Lognormal mean","Poisson Lognormal s.d.")
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
  c0 <- 1-sum(dpln(v=probv,x=0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dpln(v=probv,x=xr,approxlim=approxlim) / c0
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
llplnall <- function(v,x,cutoff=1,cutabove=1000,np=2,approxlim=10){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 aaa <- sum(nc*log(nc/n),na.rm=TRUE)+(n-sum(nc))*log((n-sum(nc))/n)+llpln(v=v,x=x,cutoff=cutoff,cutabove=cutabove,approxlim=approxlim)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
llpln <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE,
                  approxlim=10){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
 probv <- v
#probv[2] <- log(probv[2])
#if(probv[1]<=1.01 | probv[2]<=-1 | probv[1] > 50 | probv[2] > 100){
 if(probv[2]< 0.0){
  out <- -10^10
 }else{
 n <- length(x)
 if(cutoff>0){
  cprob <- 1 - sum(dpln(v=probv,x=0:(cutoff-1)))
 }else{
  cprob <- 1
 }
#
 out <- -10^10
 if(hellinger){
  tr <- 0:max(x)
  xr <- tr[tr >= cutoff]
  pdf <- dpln(v=probv,x=xr)
# pdf <- pdf / cprob
  tx <- tabulate(x+1)/length(x)
  pc <- tx[tr>=cutoff]
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  xv <- sort(unique(x))
  xp <- as.vector(table(x))
  out <- sum(xp*ldpln(v=probv,x=xv,approxlim=approxlim))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 }
 out
}
#
# Next routines to account for rounding
#
# First rounded Poisson Lognormal
#
llrpln <- function(v,x,cutoff=1,cutabove=1000,xr=1:10000,hellinger=FALSE){
 x <- x[x >= cutoff]
 x <- x[x <= cutabove]
#if(v[1] <= 2 | v[2]<= 0 | v[2] >= 1 ){
 if(v[2]< 0){
  out <- -10^6
 }else{
 probv <- v
#probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 n <- length(x)
 if(cutoff>0){
  xc <- 0:(cutoff-1)
  aaa <- dpln(v=probv,x=xc)
  cprob <- 1 - sum(aaa)
 }else{
  cprob <- 1
 }
 if(cutabove<1000){
  xc <- 0:cutabove
  aaa <- dpln(v=probv,x=xc)
  cprob <- cprob - 1 + sum(aaa)
 }
#
 out <- -10^6
  tr <- 0:max(x)
# pdf <- dpln(v=probv,x=tr[-1])
  pdf <- dpln(v=probv,x=1:800)
  pdf[ 5] <- sum(dpln(v=probv,x= 4:  6))
  pdf[10] <- sum(dpln(v=probv,x= 7: 13))
  pdf[12] <- sum(dpln(v=probv,x= 8: 16))
  pdf[15] <- sum(dpln(v=probv,x=12: 18))
  pdf[20] <- sum(dpln(v=probv,x=16: 24))
  pdf[25] <- sum(dpln(v=probv,x=20: 29))
  pdf[30] <- sum(dpln(v=probv,x=24: 36))
  pdf[35] <- sum(dpln(v=probv,x=25: 44))
  pdf[40] <- sum(dpln(v=probv,x=30: 49))
  pdf[45] <- sum(dpln(v=probv,x=35: 54))
  pdf[50] <- sum(dpln(v=probv,x=35: 64))
  pdf[55] <- sum(dpln(v=probv,x=40: 69))
  pdf[60] <- sum(dpln(v=probv,x=45: 74))
  pdf[65] <- sum(dpln(v=probv,x=50: 79))
  pdf[70] <- sum(dpln(v=probv,x=55: 84))
  pdf[75] <- sum(dpln(v=probv,x=60: 89))
  pdf[80] <- sum(dpln(v=probv,x=65: 94))
  pdf[85] <- sum(dpln(v=probv,x=70: 99))
  pdf[90] <- sum(dpln(v=probv,x=75:104))
  pdf[95] <- sum(dpln(v=probv,x=80:109))
  pdf[100] <- sum(dpln(v=probv,x=75:124))
  pdf[120] <- sum(dpln(v=probv,x=95:144))
  pdf[130] <- sum(dpln(v=probv,x=105:154))
  pdf[150] <- sum(dpln(v=probv,x=115:184))
  pdf[200] <- sum(dpln(v=probv,x=150:249))
  pdf[300] <- sum(dpln(v=probv,x=250:349))
  pdf[560] <- sum(dpln(v=probv,x=500:639))
  pdf[800] <- sum(dpln(v=probv,x=500:1100))
  pc <- tabulate(x+1, nbins=length(pdf)+1)/n
  pl <- 1:length(pdf)
  pc  <-  pc[-1][pl >= cutoff & pl <= cutabove]
  pdf <- pdf[pl >= cutoff & pl <= cutabove]
  pc[is.na(pc)] <- 0
 if(hellinger){
  pc <- pc / sum(pc,na.rm=TRUE)
  out <- -sum(((sqrt(pc) - sqrt(pdf))^2)[pc > 0])
  out <- out - 0.5 * (1-sum(pdf[pc>0])) # penalized Hellinger Distance
 }else{
  out <- n*sum(pc*log(pdf))-n*log(cprob)
 }
 if(is.infinite(out)){out <- -10^6}
 if(is.na(out)){out <- -10^6}
 }
 out
}
#
# Complete data rounded log-likelihoods
#
llrplnall <- function(v,x,cutoff=1,cutabove=1000,np=2){
 x <- x[x<=cutabove]
 n <- length(x)
 tx <- tabulate(x+1)
 tr <- 0:max(x)
 names(tx) <- paste(tr)
 nc <- tx[tr<cutoff]
 nc <- nc[nc>0]
 aaa <- sum(nc*log(nc/n))+(n-sum(nc))*log((n-sum(nc))/n)+llrpln(v=v,x=x,cutoff=cutoff,cutabove=cutabove)
 np <- np + cutoff
 aaa <- c(np,aaa,-2*aaa+np*2+2*np*(np+1)/(n-np-1),-2*aaa+np*log(n))
 names(aaa) <- c("np","log-lik","AICC","BIC")
 aaa
}
#
# Calculate the rounded Poisson Lognormal law MLE
#
rplnmle <-function(x,cutoff=1,cutabove=1000,guess=c(-1,1),
                   method="BFGS", conc=FALSE, hellinger=FALSE, hessian=TRUE){
 if(missing(guess) & hellinger){
  guess <- rplnmle(x=x,cutoff=cutoff,cutabove=cutabove,
   method=method, guess=guess,
   conc=FALSE,hellinger=FALSE)$theta
 }
 if(sum(x>=cutoff & x <= cutabove) > 0){
  aaa <- optim(par=guess,fn=llrpln,
   method=method,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  aaanm <- optim(par=guess,fn=llrpln,
   hessian=hessian,control=list(fnscale=-10),
   x=x,cutoff=cutoff,cutabove=cutabove, hellinger=hellinger)
  if(aaanm$value > aaa$value){aaa<-aaanm}
  names(aaa$par) <- c("Poisson Lognormal mean","Poisson Lognormal s.d.")
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
#probv[2] <- (probv[1]-2)/probv[2] - (probv[1]-1)
 xr <- 1:10000
 xr <- xr[xr >= cutoff]
 if(cutoff>0){
  c0 <- 1-sum(dpln(v=probv,x=0:(cutoff-1)))
 }else{
  c0 <- 1
 }
 pdf <- dpln(v=probv,x=xr) / c0
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
# Bootstrap CI for rounded Poisson Lognormal
#
bootstraprpln <- function(x,cutoff=1,cutabove=1000,
                          m=200,alpha=0.95,guess=c(-1,1),hellinger=FALSE,
                          mle.meth="rplnmle"){
#if(mle.meth!="rplnmlef"){
 aaa <- rplnmle(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}else{
# aaa <- rplnmlef(x=x,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
#}
bmles <- matrix(0,nrow=m,ncol=2)
for(i in seq(along=bmles[,1])){
 xsamp <- sample(x,size=length(x),replace=TRUE)
# if(mle.meth!="rplnmlef"){
  bmles[i,] <- rplnmle(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
# }else{
#  bmles[i,] <- rplnmlef(x=xsamp,cutoff=cutoff,cutabove=cutabove,guess=guess,hellinger=hellinger)$theta
# }
}
#
rbind(
c(quantile(bmles[,1],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa),
c(quantile(bmles[,2],c(0.5*(1-alpha),0.5,0.5+0.5*alpha),na.rm=TRUE),aaa)
)
}
#
# Compute the Poisson-Lognormal PMF
#
dpln1 <- function(v,x,cutoff=0,points=25,approxlim=10){
  quad <- gauss.hermite(points,iter=1000)
  y <- x
  xr <- x[x<=approxlim]
  gha <- exp(-exp(outer(quad[,1]*v[2]+v[1],xr*v[2]*v[2],"+")))
  y[x<=approxlim] <- exp(xr*v[1]+0.5*xr*xr*v[2]*v[2])*
                  ( quad[,2] %*% gha )/gamma(xr+1)
  xr <- x[x>approxlim]
  xrs <- (log(xr)-v[1])/v[2]
  xrb <- xrs*xrs + log(xr) - v[1] - 1
  y[x>approxlim] <- exp(-0.5*xrs*xrs)*(1+xrb/(2*xr*v[2]*v[2]))/
                                   (sqrt(2*pi)*v[2]*xr)
  if(cutoff>0){
   c0 <- 1-sum(dpln(v=v,x=0:(cutoff-1)))
   y <- y / c0
  }
  y
}
ldpln1 <- function(v,x,cutoff=0){
 log(dpln(v,x,cutoff=cutoff))
}
ldpln <- function (v, x, cutoff = 0, points = 25, approxlim = 20) 
{
    quad <- gauss.hermite(points, iterlim = 10)
    y <- x
    high <- x > approxlim
    if(any(!high)){
     xr <- x[!high]
     gha <- outer(quad[, 1]*v[2] + v[1], xr*v[2]*v[2], "+")
     gha <- exp(-exp(gha))
     aaa <- log(quad[, 2] %*% gha)
     y[!high] <- xr*v[1]+0.5*xr*xr*v[2]*v[2] + aaa - lgamma(xr + 1)
    }
    if(any(high)){
     xr <- x[high]
     xrs <- (log(xr) - v[1])/v[2]
     xrb <- xrs * xrs + log(xr) - v[1] - 1
     aaa <- log(1+xrb/(2*xr*v[2]*v[2]))
     y[high] <- -0.5*xrs*xrs + aaa - log(sqrt(2*pi)*v[2]*xr)
    }
    if (cutoff > 0) {
      c0 <- 1 - sum(dpln(v = v, x = (0:(cutoff - 1))))
      y <- y - log(c0)
    }
    y
}
dpln.refined <- function (v, x, cutoff = 0, points = 25, approxlim = 10) 
{
    quad <- gauss.hermite(points, iter = 1000)
    y <- x
    high <- x > approxlim
    if(any(!high)){
     xr <- x[!high]
     yr <- y[!high]
     for(i in seq(along=xr)){
      xri <- xr[i]
      gha <- exp(-exp(quad[,1]*v[2]+v[1]+xri*v[2]*v[2]))
      yr[i] <- exp(xri*v[1]+xri*xri*v[2]*v[2]*0.5)*sum(gha*quad[,2])/gamma(xri+1)
     }
     y[!high] <- yr
    }
    if(any(high)){
     xr <- x[high]
     xrs <- (log(xr) - v[1])/v[2]
     xrb <- xrs * xrs + log(xr) - v[1] - 1
     aaa <- log(1+xrb/(2*xr*v[2]*v[2]))
     y[high] <- exp(-0.5*xrs*xrs + aaa - log(sqrt(2*pi)*v[2]*xr))
    }
    if (cutoff > 0) {
      c0 <- 1 - sum(dpln(v = v, x = rep(0:(cutoff - 1),5)))/5
      y <- y / c0
    }
    y
}
#prpoi <- function(r,mu,sigma,points=10){
#  quad <- gauss.hermite(points)
#  gha <- exp(-exp(quad[,1]*sigma+mu+r*sigma*sigma))
#  (r*mu+r*r*sig2*0.5)*sum(gha*quad[,2])/gamma(r+1)
#}
#ldpln <- function(v,x,cutoff=0){
# log(dpln(v,x,cutoff=cutoff))
#}
dpln <- function(v,x,cutoff=0,approxlim=20){
 exp(ldpln(v,x,cutoff=cutoff,approxlim=approxlim))
}
#
# These are the basic simulation routines
#
simpln <- function(n=100, v=c(-1,1), maxdeg=10000){
  sample(x=1:maxdeg, size=n, replace=TRUE, prob=dpln(v=v,x=1:maxdeg,cutoff=1))
}
