# Load data and clusters

library(dbscan)
library(scater)
library(ggplot2)

sce <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_all_clusters.rds")

# Select conditions to analyse
conds <- c("Early Spermatocytes 1", 
           "Early Spermatocytes 2", "Mid Spermatocytes 1" , 
           "Mid Spermatocytes 2", "Late Spermatocytes 1" ,
           "Late Spermatocytes 2", "Meiosis", "S1", "S2",
           "S3", "S4", "S5", "S6", "S7", "S8", "S9",
           "S10", "S11", "S12", "S13",
           "S14")

sce <- sce[,colData(sce)$AnnotatedClusters %in% conds]

# Remove outlying cells from datasets
# B6
sce.B6.1 <- sce[,colData(sce)$Sample == "B6_do17815"]
sce.B6.1 <- normalize(sce.B6.1)
sce.B6.1 <- runTSNE(sce.B6.1)
sce.B6.1 <- runPCA(sce.B6.1)
dbs <- dbscan(reducedDims(sce.B6.1)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.B6.1)$TSNE[,1],
                  tsne2 = reducedDims(sce.B6.1)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.B6.1 <- sce.B6.1[,dbs$cluster == 1]

sce.B6.2 <- sce[,colData(sce)$Sample == "B6_do17816"]
sce.B6.2 <- normalize(sce.B6.2)
sce.B6.2 <- runTSNE(sce.B6.2)
sce.B6.2 <- runPCA(sce.B6.2)
dbs <- dbscan(reducedDims(sce.B6.2)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.B6.2)$TSNE[,1],
                  tsne2 = reducedDims(sce.B6.2)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.B6.2 <- sce.B6.2[,dbs$cluster == 1]

# Tc0
sce.Tc0.1 <- sce[,colData(sce)$Sample == "Tc0_do15984"]
sce.Tc0.1 <- normalize(sce.Tc0.1)
sce.Tc0.1 <- runTSNE(sce.Tc0.1)
sce.Tc0.1 <- runPCA(sce.Tc0.1)
dbs <- dbscan(reducedDims(sce.Tc0.1)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc0.1)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc0.1)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc0.1 <- sce.Tc0.1[,dbs$cluster == 1]

sce.Tc0.2 <- sce[,colData(sce)$Sample == "Tc0_do17622"]
sce.Tc0.2 <- normalize(sce.Tc0.2)
sce.Tc0.2 <- runTSNE(sce.Tc0.2)
sce.Tc0.2 <- runPCA(sce.Tc0.2)
dbs <- dbscan(reducedDims(sce.Tc0.2)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc0.2)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc0.2)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc0.2 <- sce.Tc0.2[,dbs$cluster == 1]

# Tc1
sce.Tc1.1 <- sce[,colData(sce)$Sample == "Tc1_do15983"]
sce.Tc1.1 <- normalize(sce.Tc1.1)
sce.Tc1.1 <- runTSNE(sce.Tc1.1)
sce.Tc1.1 <- runPCA(sce.Tc1.1)
dbs <- dbscan(reducedDims(sce.Tc1.1)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc1.1)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc1.1)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc1.1 <- sce.Tc1.1[,dbs$cluster == 1]

sce.Tc1.2 <- sce[,colData(sce)$Sample == "Tc1_do17623"]
sce.Tc1.2 <- normalize(sce.Tc1.2)
sce.Tc1.2 <- runTSNE(sce.Tc1.2)
sce.Tc1.2 <- runPCA(sce.Tc1.2)
dbs <- dbscan(reducedDims(sce.Tc1.2)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.Tc1.2)$TSNE[,1],
                  tsne2 = reducedDims(sce.Tc1.2)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.Tc1.2 <- sce.Tc1.2[,dbs$cluster == 1]

# Juvenile
sce.P30 <- sce[,colData(sce)$Sample == "P30_do17825"]
sce.P30 <- normalize(sce.P30)
sce.P30 <- runTSNE(sce.P30)
sce.P30 <- runPCA(sce.P30)
dbs <- dbscan(reducedDims(sce.P30)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.P30)$TSNE[,1],
                  tsne2 = reducedDims(sce.P30)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.P30 <- sce.P30[,dbs$cluster == 1]

sce.P35 <- sce[,colData(sce)$Sample == "P35_do17827"]
sce.P35 <- normalize(sce.P35)
sce.P35 <- runTSNE(sce.P35)
sce.P35 <- runPCA(sce.P35)
dbs <- dbscan(reducedDims(sce.P35)$PCA, eps = 1.5)
ggplot(data.frame(tsne1 = reducedDims(sce.P35)$TSNE[,1],
                  tsne2 = reducedDims(sce.P35)$TSNE[,2],
                  clusters = as.factor(dbs$cluster))) + geom_point(aes(tsne1, tsne2, colour = clusters))
sce.P35 <- sce.P35[,dbs$cluster == 1]

# Split into matrices
B6_1 <- as.matrix(counts(sce.B6.1))
B6_2 <- as.matrix(counts(sce.B6.2))

Tc0_1 <- as.matrix(counts(sce.Tc0.1))
Tc0_2 <- as.matrix(counts(sce.Tc0.2))

Tc1_1 <- as.matrix(counts(sce.Tc1.1))
Tc1_2 <- as.matrix(counts(sce.Tc1.2))

P30 <- as.matrix(counts(sce.P30))
P35 <- as.matrix(counts(sce.P35))

# Rename samples
colnames(B6_1) <- paste("B6.1", colData(sce.B6.1)$AnnotatedClusters,
                        colData(sce.B6.1)$Barcode, sep = "_")
colnames(B6_2) <- paste("B6.2", colData(sce.B6.2)$AnnotatedClusters,
                        colData(sce.B6.2)$Barcode, sep = "_")

colnames(Tc0_1) <- paste("Tc0.1", colData(sce.Tc0.1)$AnnotatedClusters,
                        colData(sce.Tc0.1)$Barcode, sep = "_")
colnames(Tc0_2) <- paste("Tc0.2", colData(sce.Tc0.2)$AnnotatedClusters,
                        colData(sce.Tc0.2)$Barcode, sep = "_")

colnames(Tc1_1) <- paste("Tc1.1", colData(sce.Tc1.1)$AnnotatedClusters,
                        colData(sce.Tc1.1)$Barcode, sep = "_")
colnames(Tc1_2) <- paste("Tc1.2", colData(sce.Tc1.2)$AnnotatedClusters,
                        colData(sce.Tc1.2)$Barcode, sep = "_")

colnames(P30) <- paste("Juvenile.1", colData(sce.P30)$AnnotatedClusters,
                         colData(sce.P30)$Barcode, sep = "_")
colnames(P35) <- paste("Juvenile.2", colData(sce.P35)$AnnotatedClusters,
                         colData(sce.P35)$Barcode, sep = "_")


# Merge datasets
B6 <- cbind(B6_1, B6_2)
B6 <- B6[rowMeans(B6) > 0.1,]
Tc0 <- cbind(Tc0_1, Tc0_2)
Tc0 <- Tc0[rowMeans(Tc0) > 0.1,]
Tc1 <- cbind(Tc1_1, Tc1_2)
Tc1 <- Tc1[rowMeans(Tc1) > 0.1,]
Juvenile <- cbind(P30, P35)
Juvenile <- Juvenile[rowMeans(Juvenile) > 0.1,]

saveRDS(B6, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_B6.rds")
saveRDS(Tc0, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_Tc0.rds")
saveRDS(Tc1, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_Tc1.rds")
saveRDS(Juvenile, "Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/All_Juvenile.rds")
