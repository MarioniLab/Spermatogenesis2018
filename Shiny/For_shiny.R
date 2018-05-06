#### This script collects all datasets and prepaes them to be loaded into shiny

library(scater)
library(scran)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

# Read in data
# Adult
sce <- readRDS("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/10X_data/SCE_all_clusters.rds")

# Read in genanames
mouse.genes <- read.table("Dropbox (Cambridge University)/SST_spermatocytes/Analysis/data/Mouse_genes.txt",
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mouse.genes <- mouse.genes[!grepl("CHR", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("GL", mouse.genes$Chromosome.scaffold.name) &
                             !grepl("JH", mouse.genes$Chromosome.scaffold.name),]
mouse.genes$Chromosome.scaffold.name <- paste("Chr", 
                                              mouse.genes$Chromosome.scaffold.name,
                                              sep = "")

# Create colour vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_vector = c("#242435", "#D9D9ED","#D62261","#BCBCC2","#B19B31","#282B69","#572972","#421449","#4A1028","#851A4F","#643F18","#0F1131","#FFDFB2","#60AC78","#985D25","#176D38","#C58C58","#827BAF","#F36E21","#71719A","#47476B","#AFB3DA","#9595C9","#93381F","#A58D81","#8A8A9B","#2A348B","#6167AF","#5F4E97","#EE6493","#CECAE5","#ff00ff")
names(color_vector) <- unique(colData(sce)$Cluster)

# Save output
save.image("Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce.RData")
