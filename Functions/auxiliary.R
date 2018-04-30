###############################################
#### Script to collect auxiliary functions ####
###############################################

# Load libraries needed in this script
library(scater)
library(scran)
library(dynamicTreeCut)
library(slingshot)
library(dbscan)

#### Split Single cell experiment
split.sce <- function(sce, groups, colData.name = "Cluster"){
  # List to collect individual single cell experiments
  list.out <- list()
  for(i in groups){
    cur_sce <- sce[,colData(sce)[[colData.name]] == i]
    cur_sce <- normalize(cur_sce)
    list.out[[i]] <- cur_sce
  }
  names(list.out) <- groups
  list.out
}


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
  else if(typeof(sce) == "list") {
    lapply(sce, function(n){
      HVG <- trendVar(n, use.spikes = FALSE)
      HVG.1 <- decomposeVar(n, HVG)
      HVG.1 <- HVG.1[order(HVG.1$bio, decreasing = TRUE),]
      rownames(HVG.1)[1:numberGenes]
    })
  }
  else{
    print("First argument must either be a single cell experiment object \n
          or a list")
  }
}

#### Clustering
# Perform clustering using dynamic tree cut
DTC <- function(sce, HVG.genes, minClusterSize = 10, deepSplit = 0){
  if(typeof(sce) == "S4"){
    dist.all <- as.dist(sqrt((1 - cor(as.matrix(logcounts(sce)[HVG.genes,]), 
                                      method = "spearman"))/2))
    
    dendro <- hclust(dist.all, method = "ward.D2")
    
    ct <- as.character(cutreeDynamic(dendro = dendro, distM = as.matrix(dist.all), 
                        minClusterSize = minClusterSize, deepSplit = deepSplit))
  }
  else {
    out <- list()
    for(i in 1:length(sce)){
      dist.all <- as.dist(sqrt((1 - cor(as.matrix(logcounts(sce[[i]])[HVG.genes[[i]],]), 
                                        method = "spearman"))/2))
      
      dendro <- hclust(dist.all, method = "ward.D2")
      
      cur_clusters <- paste(names(sce)[i], as.character(cutreeDynamic(dendro = dendro, 
                                                      distM = as.matrix(dist.all), 
                minClusterSize = minClusterSize, deepSplit = deepSplit)), sep = "_")
      names(cur_clusters) <- colData(sce[[i]])$Barcode
      
      out[[names(sce)[i]]] <- cur_clusters
    }
    names(out) <- names(sce)
    out
  }
}

#### Find specifc marker genes
marker.detection <- function(sce, clusters){
  # User scran function findMarkers to perform differential expression
  cur_markers <- findMarkers(sce, clusters)
  
  # Collect group specific markers
  markers.spec <- lapply(cur_markers, function(n){
    cur_n <- n[n$FDR < 0.1 & apply(n[,3:ncol(n)], 1, function(x){sum(x > 0)}) == ncol(n) - 2,]
    cur_n$GeneName <- rowData(sce)$Symbol[match(rownames(cur_n), rowData(sce)$ID)]
    cur_n
  })
  
}

#### Compute pseudorank
PT <- function(rd, clusters, col_vector, 
               exclude = NULL, start = NULL, end = NULL){
    if(!is.null(exclude)){
      cur_rd <- rd[!exclude,]
      
      cur_lin <- principal.curve(cur_rd)
      
      plot(cur_rd, col = col_vector[clusters[!exclude]], 
           pch = 16, type = "p")
      lines(cur_lin, lwd = 3)
      
      mat.out <- matrix(data = NA, ncol = ncol(cur_rd) + 1, nrow = length(clusters))
      rownames(mat.out) <- names(clusters)
      colnames(mat.out) <- c(colnames(cur_rd), "rank")
      
      mat.out[!exclude,1:ncol(cur_rd)] <- cur_lin$s
      mat.out[!exclude,"rank"] <- order(cur_lin$tag)
      
      mat.out
    }
    else{
      cur_rd <- rd
      
      cur_lin <- principal.curve(cur_rd)
      
      plot(cur_rd, col = col_vector[clusters], 
           pch = 16, type = "p")
      lines(cur_lin, lwd = 3)
      
      mat.out <- matrix(data = NA, ncol = ncol(cur_rd) + 1, nrow = length(clusters))
      rownames(mat.out) <- names(clusters)
      colnames(mat.out) <- c(colnames(cur_rd), "rank")
      
      mat.out[!exclude,1:ncol(cur_rd)] <- cur_lin$s
      mat.out[!exclude,"rank"] <- order(cur_lin$tag)
      
      mat.out
    }
}
