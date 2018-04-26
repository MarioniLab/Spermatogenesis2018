###############################################
#### Script to collect auxiliary functions ####
###############################################

#### HVG
# Calculate highly variable genes
# sce: one or multiple single cell experiment objects in form of a list
HVG <- function(sce, numberGenes = 1000){
  # One single cell expereiment object
  if(typeof(sce) == "S4"){
    HVG <- trendVar(sce, use.spikes = FALSE)
    HVG.1 <- decomposeVar(sce, HVG)
    HVG.1 <- HVG.1[order(HVG.1$bio, decreasing = TRUE),]
    rownames(HVG.1)[1:numberGenes]
  }
  # Multiple single cell experiment objects
  else {
    lapply(sce, function(n){
      HVG <- trendVar(n, use.spikes = FALSE)
      HVG.1 <- decomposeVar(n, HVG)
      HVG.1 <- HVG.1[order(HVG.1$bio, decreasing = TRUE),]
      rownames(HVG.1)[1:numberGenes]
    })
  }
}

#### Clustering
# Perform clustering using dynamic tree cut
DTC <- function(sce, HVG.genes, minClusterSize = 10, deepSplit = 0){
  if(typeof(sce) == "S4"){
    dist.all <- as.dist(sqrt((1 - cor(as.matrix(logcounts(sce)[HVG.genes,]), 
                                      method = "spearman"))/2))
    
    dendro <- hclust(dist.all, method = "ward.D2")
    
    ct <- cutreeDynamic(dendro = dendro, distM = as.matrix(dist.all), 
                        minClusterSize = minClusterSize, deepSplit = deepSplit)
  }
  else {
    out <- list()
    for(i in 1:length(sce)){
      dist.all <- as.dist(sqrt((1 - cor(as.matrix(logcounts(sce[[i]])[HVG.genes[[i]],]), 
                                        method = "spearman"))/2))
      
      dendro <- hclust(dist.all, method = "ward.D2")
      
      out[[names(sce)[i]]] <- cutreeDynamic(dendro = dendro, 
                          distM = as.matrix(dist.all), 
                          minClusterSize = minClusterSize, deepSplit = deepSplit)
    }
  }
}
