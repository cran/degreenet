#
# Find a graph from a given degree sequence
# but not a random graph
#
reedmolloy <- function(deg,
                       greedy=FALSE,
                       verbose=TRUE){
 n <- length(deg) 
 mdeg <- max(deg) 
 while(2*floor(sum(deg)/2) != sum(deg) ){
  if(verbose & sum(deg)>1){
   stop(paste("You have provided an odd number of links =",sum(deg),"\n",
       "This does not make sense for this network, and so resample...\n"))
  }
  if(verbose){print(table(deg))}
  mdeg <- max(deg)
 }
 deg <- deg[order(-deg)]
 sm <- matrix(0,ncol=n,nrow=n)
 i <- 1
 while(i <= n){
# if(verbose){print(i)}
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
  stop("need repeat sample")
 }
 if(require(network, quietly=TRUE)){
  network(sm, directed=FALSE)
 }else{
  sm
 }
}
