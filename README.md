# Transcriptional dynamics during juvenile and adult spermatogenesis

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq, bulk RNAseq and CUT&RUN data.

## Obtaining the data

The analyses scripts load in a `SingleCellExperiment` object that contains either the cell ranger filtered cells or the EmptyDrops filtered cells (see Methods section of the publication).
These single-cell RNA sequencing data have been deposited on ArrayExpress under the accession number [E-MTAB-6946](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6946/).
The following code chunk downloads the data and stores them in a `SingleCellExperiment` object that can be used to reproduce the analysis performed for this project.

In the first step, clone this repository to your local Github folder:

```{bash}
cd ~/Github
git clone https://github.com/MarioniLab/Spermatogenesis2018.git
```

To obtain the cell ranger filtered data use:

```{r}
setwd("~/GitHub/Spermatogenesis2018/")

# Raw cellranger filtered transcriptomes
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/E-MTAB-6946.processed.1.zip", 
                destfile = "cellranger_raw.zip")
unzip("cellranger_raw.zip")

# Cellranger filtered metadata
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/E-MTAB-6946.processed.2.zip", 
                destfile = "cellranger_metadata.zip")
unzip("cellranger_metadata.zip")  
                
# Create SingleCellExperiment object
library(SingleCellExperiment)
library(Matrix)

cellranger_raw <- readMM("raw_counts.mtx")
cellranger_metadata <- read.table("cell_metadata.txt", sep = " ")
sce_cellranger <- SingleCellExperiment()  
```

To obtain the EmptyDrops filtered data use:

```{r}
# Raw empty drops filtered transcriptomes
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/E-MTAB-6946.processed.3.zip", 
                destfile = "EmptyDrops_raw.zip")
unzip("EmptyDrops_raw.zip")
               
# EmptyDrops filtered metadata
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/E-MTAB-6946.processed.4.zip", 
                destfile = "EmptyDrops_metadata.zip")
unzip("EmptyDrops_metadata.zip")                

sce_emptyDrops <- 
```

## Analysis scripts

The following folders contain scripts for data processing and analysis. 
A short description can be found below:

* [Preprocessing](../master/Preprocessing/) contains quality control and normalization scripts for scRNA-Seq and bulk RNA-Seq data, as well as a scropt to process the data from Chen et al., Cell Research, 2018.
 
* All scripts to reproduce the individual figures (main and supplements) can be found in [Figures](../master/Figures).

* [Data](../master/Data) contains a file listing mouse gene IDs, Symbols and corresponding chromosome information.

* [Functions](../master/Functions) contains an R script with wrapper functions for simple tasks throughout the analysis. Furthermore, this folder contains a script to count and normalize reads of CUT&RUN and ChIP-Seq libraries in promoter regions.

* [Analysis](../master/Analysis) contains .Rmd files with additional analysis performed on the data. These analyses have not been used to present data in the main manuscript.

* [Shiny](../master/Shiny) conatins scripts to prepare data for the shiny app visualization and the server.R and ui.R files.

* [Shiny_ED](../master/Shiny_ED): Same as Shiny but for emptyDrops data .

The preprint can be found at:

[Staged developmental mapping and X chromosome transcriptional dynamics during mouse spermatogenesis](https://www.biorxiv.org/content/early/2018/06/20/350868)
