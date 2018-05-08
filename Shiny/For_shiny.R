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
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_vector = c("#47476B", "#F36E21","#985D25","#643F18","#5F4E97","#C58C58","#B19B31","#572972","#D62261","#851A4F","#FFDFB2","#71719A","#9595C9","#AFB3DA","#176D38","#421449","#EE6493","#2A348B","#827BAF","#93381F","#D9D9ED","#242435","#CECAE5","#067277","#8A8A9B","#6167AF","#0F1131","#282B69","#BCBCC2","#4A1028","#A58D81","#ff00ff")
names(color_vector) <- levels(as.factor(colData(sce)$Clusters))

# Save output
save.image("Dropbox (Cambridge University)/SST_spermatocytes/Shiny/data/sce.RData")
