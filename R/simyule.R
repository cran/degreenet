#
# Find a graph from a Yule distribution
# but not a random graph
#
ryule <- function(n=20,rho=2.5,cutoff=1,cutabove=1000,
                  greedy=FALSE,
                  maxdeg=10000,maxout=TRUE,verbose=FALSE){
 mdeg <- n 
 deg <- 1
 while(mdeg>n-1 | 2*floor(sum(deg)/2) != sum(deg) ){
  if(verbose & sum(deg)>1){
   print(paste("odd number of links =",sum(deg),", so resampling..."))}
  deg <- sample(x=1:maxdeg,size=n,replace=TRUE,prob=dyule(v=rho,x=1:maxdeg))
  mdeg <- max(deg)
  if(verbose & mdeg>n-1){print(paste("max.deg =",mdeg,"> n - 1, so resampling..."))}
  if(verbose){print(table(deg))}
  if(maxout){deg[deg>n-1] <- n-1}
  mdeg <- max(deg)
 }
 deg <- deg[order(-deg)]
 sm <- matrix(0,ncol=n,nrow=n)
 i <- 1
 while(i <= n){
  if(deg[i] <= sum(deg[-i]>0)){
   if(deg[i]>0){
    if(greedy){
     x <- (((1:n)[-i])[deg[-i]>0])[1:deg[i]]
    }else{
     x <- sample(x=(1:n)[-i],size=deg[i],prob=deg[-i]/sum(deg[-i]))
    }
    deg[x] <- deg[x] - 1
    deg[i] <- 0
    sm[i,x] <- 1
    sm[x,i] <- 1
   }
  }else{
   if(verbose){
     print(paste("no fit on",i,"th node with",sum(deg),"links unallocated"))}
   if(verbose){print(table(deg))}
   i <- n + 2
  }
  i <- i + 1
 }
 if(sum(deg)>0 | i > n + 1){
  warning("repeat sample")
  sm <- ryule(n=n,rho=rho,cutoff=cutoff,cutabove=cutabove,maxdeg=maxdeg,maxout=maxout)
 }
 if(require(network, quietly=TRUE)){
  network(sm, directed=FALSE)
 }else{
  sm
 }
}
