# Load data and clusters

library(dbscan)
library(scater)
library(ggplot2)

sce <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_all.rds")

# Select conditions to analyse
sce <- sce[,colData(sce)$AnnotatedClusters %in% levels(colData(sce)$AnnotatedClusters)[1:22]]

# Remove outlying cells from datasets
# B6
sce.B6 <- sce[,colData(sce)$Sample == "B6"]
sce.B6 <- normalize(sce.B6)
sce.B6 <- runPCA(sce.B6)
dbs <- dbscan(reducedDims(sce.B6)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.B6)$TSNE[,1],
                  tsne2 = reducedDims(sce.B6)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.B6 <- sce.B6[,dbs$cluster == 1]

# Tc0
sce.Tc0 <- sce[,colData(sce)$Sample == "Tc0"]
sce.Tc0 <- normalize(sce.Tc0)
sce.Tc0 <- runPCA(sce.Tc0)
dbs <- dbscan(reducedDims(sce.Tc0)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc0)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc0)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc0 <- sce.Tc0[,dbs$cluster == 1]

# Tc1
sce.Tc1 <- sce[,colData(sce)$Sample == "Tc1"]
sce.Tc1 <- normalize(sce.Tc1)
sce.Tc1 <- runPCA(sce.Tc1)
dbs <- dbscan(reducedDims(sce.Tc1)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc1)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc1)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc1 <- sce.Tc1[,dbs$cluster == 1]

# Split into matrices
B6 <- as.matrix(counts(sce.B6))

Tc0 <- as.matrix(counts(sce.Tc0))

Tc1 <- as.matrix(counts(sce.Tc1))

# Rename samples
colnames(B6) <- paste(paste(colData(sce.B6)$Sample, colData(sce.B6)$Library, sep = ""), 
                        colData(sce.B6)$AnnotatedClusters,
                        colData(sce.B6)$Barcode, sep = "-")

colnames(Tc0) <- paste(paste(colData(sce.Tc0)$Sample, colData(sce.Tc0)$Library, sep = ""), 
                      colData(sce.Tc0)$AnnotatedClusters,
                      colData(sce.Tc0)$Barcode, sep = "-")

colnames(Tc1) <- paste(paste(colData(sce.Tc1)$Sample, colData(sce.Tc1)$Library, sep = ""), 
                      colData(sce.Tc1)$AnnotatedClusters,
                      colData(sce.Tc1)$Barcode, sep = "-")


# Merge datasets
B6.final <- B6[rowMeans(B6) > 0.1 & rowMeans(Tc0) > 0.1 & rowMeans(Tc1) > 0.1,]
Tc0.final <- Tc0[rowMeans(B6) > 0.1 & rowMeans(Tc0) > 0.1 & rowMeans(Tc1) > 0.1,]
Tc1.final <- Tc1[rowMeans(B6) > 0.1 & rowMeans(Tc0) > 0.1 & rowMeans(Tc1) > 0.1,]


saveRDS(B6.final, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_B6.rds")
saveRDS(Tc0.final, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_Tc0.rds")
saveRDS(Tc1.final, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_Tc1.rds")
