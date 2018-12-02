# Transcriptional dynamics during juvenile and adult spermatogenesis

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq, bulk RNAseq and CUT&RUN data.

* [Preprocessing](../master/Preprocessing/) contains quality control and normalization scripts for scRNA-Seq and bulk RNA-Seq data, as well as a scropt to process the data from Chen et al., Cell Research, 2018.
 
* All scripts to reproduce the individual figures (main and supplements) can be found in [Figures](../master/Figures).

* [Data](../master/Data) contains a file listing mouse gene IDs, Symbols and corresponding chromosome information.

* [Functions](../master/Functions) contains an R script with wrapper functions for simple tasks throughout the analysis. Furthermore, this folder contains a script to count and normalize reads of CUT&RUN and ChIP-Seq libraries in promoter regions.

* [Analysis](../master/Analysis) contains .Rmd files with additional analysis performed on the data. These analyses have not been used to present data in the main manuscript.

* [Shiny](../master/Shiny) conatins scripts to prepare data for the shiny app visualization and the server.R and ui.R files.

* [Shiny_ED](../master/Shiny_ED): Same as Shiny but for emptyDrops data .

The preprint can be found at:

[Staged developmental mapping and X chromosome transcriptional dynamics during mouse spermatogenesis](https://www.biorxiv.org/content/early/2018/06/20/350868)
