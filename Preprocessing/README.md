# Scripts to process the data for analysis

To recreate the preprocessing done in the paper run the scripts in following order:

* Basic quality control and filtering script located in [QualityControl](../../master/Preprocessing/QualityControl/).

* The clustering script to detect cell types in B6 located in [Clustering](../../master/Preprocessing/Clustering/) ([CellTypesB6.Rmd](../../master/Preprocessing/Clustering/CellTypesB6.Rmd))

* Scripts that perform subclustering on enriched cell populations (P10, P20 and P35) located in [Clustering](../../master/Preprocessing/Clustering/) ([P10Enrichment.Rmd](../../master/Preprocessing/Clustering/P10Enrichment.Rmd), [P20Enrichment.Rmd](../../master/Preprocessing/Clustering/P20Enrichment.Rmd)  and [P35Enrichment.Rmd](../../master/Preprocessing/Clustering/P35Enrichment.Rmd))

