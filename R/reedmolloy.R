#  File degreenet/R/reedmolloy.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
#
# Find a graph from a given degree sequence
# but not a random graph
#
reedmolloy <- function(deg,
                       greedy=FALSE,
                       warn=TRUE,
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
#sm <- matrix(0,ncol=n,nrow=n)
 sm <- NULL
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
#   sm[i,x] <- 1
#   sm[x,i] <- 1
    sm <- rbind(sm,c(max(i,x),min(i,x)))
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
  if(!warn){
    stop("need repeat sample")
  }else{
    warning("need repeat sample")
  }
 }
 if(require(network, quietly=TRUE)){
  smn <- network.initialize(n, directed=FALSE)
  smn <- add.edges(smn,as.list(sm[,1]),as.list(sm[,2]))
# network(sm, directed=FALSE)
 }else{
  sm
 }
}
