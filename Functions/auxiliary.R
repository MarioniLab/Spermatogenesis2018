###############################################
#### Script to collect auxiliary functions ####
###############################################

# Calculate highly variable genes
# sce: one or multiple single cell experiment objects in form of a list
HVG <- function(sce, numberGenes = 1000){
  # One single cell expereiment object
  if(length(sce) == 1){
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
