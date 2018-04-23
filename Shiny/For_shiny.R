#### This script collects all datasets and prepaes them to be loaded into shiny

library(scater)
library(scran)
library(Rtsne)
library(ggplot2)

# Read in data
# Adult
sce.adult <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_all.rds")
sce.adult <- sce.adult[,grepl("B6", colData(sce.adult)$Sample)]
sce.adult <- normalize(sce.adult)
rowData(sce.adult) <- rowData(sce.adult)[,1:2]
colData(sce.adult) <- colData(sce.adult)[,c("Sample", "Barcode", "Cluster")]

# Juvenile
sce.juvenile <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_juvenile.rds")
rowData(sce.juvenile) <- rowData(sce.juvenile)[,1:2]
colData(sce.juvenile) <- colData(sce.juvenile)[,c("Sample", "Barcode", "Cluster")]

# P10
sce.P10 <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_P10.rds")
rowData(sce.P10) <- rowData(sce.P10)[,1:2]
colData(sce.P10) <- colData(sce.P10)[,c("Sample", "Barcode", "Clusters")]
colData(sce.P10)$Cluster <- colData(sce.P10)$Clusters
colData(sce.P10)$Clusters <- NULL

# P35
sce.P35 <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_P35.rds")
rowData(sce.P35) <- rowData(sce.P35)[,1:2]
colData(sce.P35) <- colData(sce.P35)[,c("Sample", "Barcode", "Cluster")]

# Merge datasets
genes <- intersect(intersect(intersect(rownames(sce.adult), rownames(sce.juvenile)), 
                             rownames(sce.P10)), rownames(sce.P35))

sce <- cbind(sce.adult[genes,], sce.juvenile[genes,], sce.P10[genes,], sce.P35[genes,])
sce <- normalize(sce)

#### Remove batch effects
# Calculate highly variable genes
B6.1 <- sce[,grepl("B6_do17815", colData(sce)$Sample)]
B6.1 <- normalize(B6.1)
HVG.B6.1 <- trendVar(B6.1, use.spikes = FALSE)
HVG.B6.1 <- decomposeVar(B6.1, HVG.B6.1)

B6.2 <- sce[,grepl("B6_do17816", colData(sce)$Sample)]
B6.2 <- normalize(B6.2)
HVG.B6.2 <- trendVar(B6.2, use.spikes = FALSE)
HVG.B6.2 <- decomposeVar(B6.2, HVG.B6.2)

P10 <- sce[,grepl("P10_do17821", colData(sce)$Sample)]
P10 <- normalize(P10)
HVG.P10 <- trendVar(P10, use.spikes = FALSE)
HVG.P10 <- decomposeVar(P10, HVG.P10)

P15 <- sce[,grepl("P15_do17828", colData(sce)$Sample)]
P15 <- normalize(P15)
HVG.P15 <- trendVar(P15, use.spikes = FALSE)
HVG.P15 <- decomposeVar(P15, HVG.P15)

P20 <- sce[,grepl("P20_do17824", colData(sce)$Sample)]
P20 <- normalize(P20)
HVG.P20 <- trendVar(P20, use.spikes = FALSE)
HVG.P20 <- decomposeVar(P20, HVG.P20)



P30 <- sce[,grepl("P30_do17825", colData(sce)$Sample)]
P30 <- normalize(P30)
HVG.P30 <- trendVar(P30, use.spikes = FALSE)
HVG.P30 <- decomposeVar(P30, HVG.P30)

P35 <- sce[,grepl("P35_do17827", colData(sce)$Sample)]
P35 <- normalize(P35)
HVG.P35 <- trendVar(P35, use.spikes = FALSE)
HVG.P35 <- decomposeVar(P35, HVG.P35)

HVG.df <- combineVar(HVG.B6.1, HVG.B6.2, HVG.P10,
                     HVG.P15, HVG.P20, HVG.P30,
                     HVG.P35)

HVG.df <- HVG.df[order(HVG.df$bio, decreasing = TRUE),]
HVG <- rownames(HVG.df)[1:1000]

# Batch correction
corrected <- mnnCorrect(as.matrix(logcounts(B6.1)[HVG,]),
                        as.matrix(logcounts(B6.2)[HVG,]),
                        as.matrix(logcounts(P10)[HVG,]),
                        as.matrix(logcounts(P15)[HVG,]),
                        as.matrix(logcounts(P20)[HVG,]),
                        as.matrix(logcounts(P30)[HVG,]),
                        as.matrix(logcounts(P35)[HVG,]),
                        cos.norm.in = TRUE, cos.norm.out = TRUE,
                                  sigma = 0.1)

# Compute tSNE
set.seed(123)
tsne <- Rtsne(t(cbind(corrected$corrected[[1]], 
                      corrected$corrected[[2]],
                      corrected$corrected[[3]],
                      corrected$corrected[[4]],
                      corrected$corrected[[5]],
                      corrected$corrected[[6]],
                      corrected$corrected[[7]])))

sce <- cbind(B6.1, B6.2, P10, P15, P20, P30, P35)
sce <- normalize(sce)

# Plot tSNE
ggplot(data = data.frame(
  tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], 
  sample = factor(colData(sce)$Sample))) + 
  geom_point(aes(tsne1, tsne2, colour = sample)) + theme_minimal()


reducedDims(sce)$TSNE <- tsne$Y

# Read in genanames
mouse.genes <- read.table("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/Mouse_genes.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mouse.genes <- mouse.genes[!grepl("CHR", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("GL", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("JH", mouse.genes$Chromosome.scaffold.name),]
mouse.genes$Chromosome.scaffold.name <- paste("Chr", 
                                              mouse.genes$Chromosome.scaffold.name,
                                              sep = "")

# Save output
save.image("Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce.RData")
