#
# These are the basic simulation routines
#
simyule <- function(n=100, rho=4, maxdeg=10000){
  sample(x=1:maxdeg, size=n, replace=TRUE, prob=dyule(v=rho,x=1:maxdeg))
}
simnb <- function(n=100, v=c(5,0.2), maxdeg=10000){
  sample(x=1:maxdeg, size=n, replace=TRUE,
         prob=dnbinom(x=1:maxdeg, size=v[2]*v[1], prob=v[2]))
}
simwar <- function(n=100, v=c(3.5,0.1), maxdeg=10000){
  sample(x=1:maxdeg, size=n, replace=TRUE,
         prob=dwar(v=v,x=1:maxdeg))
}
simdp <- function(n=100, v=3.5, maxdeg=10000){
  sample(x=1:maxdeg, size=n, replace=TRUE,
         prob=ddp(v=v,x=1:maxdeg))
}
