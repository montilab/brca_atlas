#### To use set clustFunc=cKmeansWrapperSqrt in K2preproc()
library(conclust)

## Wrapper to subsample
cKmeansWrapperSqrt <- function(dataMatrix, clustList) {
  
  clustList$labs <- as.character(clustList$labs)
  
  eMatSub <- clustList$eMat[rownames(dataMatrix), clustList$labs %in%
                              colnames(dataMatrix)]
  labsSub <- clustList$labs[clustList$labs %in% colnames(dataMatrix)]
  
  ## Subsample the data
  sVec <- unlist(lapply(unique(labsSub), function(x) {
    wsamps <- which(labsSub == x)
    sample(wsamps, sqrt(length(wsamps)))
  }))
  eMatSub <- eMatSub[, sVec]
  labsSub <- labsSub[sVec]
  
  ## Get constraints
  mustLink <- outer(labsSub, labsSub, "==")
  mustLink[upper.tri(mustLink, diag=TRUE)] <- FALSE
  mustLink <- which(mustLink, arr.ind=TRUE)
  
  ## Cluster data
  dClust=factor(conclust::lcvqe(t(eMatSub), k=2, mustLink=mustLink,
                      cantLink=matrix(c(1, 1), nrow=1), maxIter=clustList$maxIter),
                levels=c(1, 2))
  
  ## Get label-level clusters
  dMat <- as.matrix(table(dClust, labsSub))[, colnames(dataMatrix)]
  modVec <- apply(dMat, 2, which.max)
  mods <- paste(modVec, collapse="")
  
  return(mods)
}