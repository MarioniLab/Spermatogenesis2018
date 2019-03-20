#### This script collects all datasets and prepaes them to be loaded into shiny

library(scater)
library(scran)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# Read in data
# Adult
sce <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_emptyDrops.rds")

# Read in genanames
mouse.genes <- read.table("GitHub/Spermatogenesis2018/Data/Mouse_genes.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mouse.genes <- mouse.genes[!grepl("CHR", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("GL", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("JH", mouse.genes$Chromosome.scaffold.name),]
mouse.genes$Chromosome.scaffold.name <- paste("Chr", 
                                              mouse.genes$Chromosome.scaffold.name,
                                              sep = "")

# Create colour vector
color_vector = metadata(sce)$color_vector

# Downsample sce object
#set.seed(12345)
#sce <- sce[,sample(1:ncol(sce), 20000)]

# Save output
save.image("Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce_ED.RData")
