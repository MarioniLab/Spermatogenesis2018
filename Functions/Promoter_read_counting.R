# This script counts reads in promoters from CUT & RUN and ChIP-Seq data
# The reads are normalized based on the total library size (CPM)

# Read in libraries
library(csaw)
library(EnsDb.Mmusculus.v79)
library(rtracklayer)

# Generate necessary annotations
prom <- promoters(genes(EnsDb.Mmusculus.v79))
prom <- prom[seqnames(prom) %in% c(as.character(1:19), "X", "Y", "MT"),]
seqlevels(prom) <- c(as.character(1:19), "X", "Y", "MT")

# Blacklisted regions
black <- import("http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz")
seqlevels(black) <- sub("chr", "", seqlevels(black))
seqnames(black) <- sub("chr", "", seqnames(black))

# Set parameters
# Single of paired-end 
reads <- "single"

# Input directory
input <- "Dropbox (Cambridge University)/SST_spermatocytes/Revisions/Data/Hammoud_CSC/ChipSeq/bam/"

# Output directory
output <- "Dropbox (Cambridge University)/SST_spermatocytes/Revisions/Data/Hammoud_CSC/ChipSeq/promoters/"


# Parameters for parsing the bam file
if(reads == "single"){
  param <- readParam(dedup=TRUE, minq=10, discard=black, pe="none")
} else{
  param <- readParam(dedup=TRUE, minq=10, discard=black, pe="both", max.frag=1000)
}

# Loop through the files, count read in promoters, normalize and save as .txt file
files <- list.files(input, full.names = TRUE, pattern = "bam$")
names <- as.character(sapply(list.files(input, full.names = FALSE, pattern = "bam$"), 
                             function(n){unlist(strsplit(n, "\\."))[1]}))

for(i in 1:length(files)){
  cur_bam <- files[i]
  prom.counts <- regionCounts(cur_bam, regions = prom, param=param)
  cur_bins <- windowCounts(cur_bam, bin = TRUE, width = 1000, param=param)
  cur_counts <- assays(prom.counts)$counts
  cur_rpm <- t(t(cur_counts)/(colSums(assays(cur_bins)$counts)/1000000))
  colnames(cur_rpm) <- names[i]
  
  write.table(cur_rpm, paste(output, names[i], ".txt", sep = ""))
}
