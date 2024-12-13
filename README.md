# diapause-sc-multiome

#### This project contains analysis of single cell multiomic data from multiple timepoints in killifish development.

TODOs:

-   High Priority
    -   [x] adjust cell calling parameters in cellranger
    -   [ ] create UMAPs for each sample and combinations of samples
        -   create plot saving function
    -   [ ] Cell annotations
        -   manual
        -   automated
            -   <https://doi.org/10.1016/j.csbj.2021.01.015>
    -   [ ] Marker gene/DE: try different methods
        -   MAST
-   Medium Priority
    -   compare with bulk data
    -   density UMAP, heatmaps
    -   ATAC UMAP
    -   Doublet identification with DF
    -   Peak-gene linkages
    -   convert ncbi gene names to param's annotations
    -   split 1_preprocessing into multiple notebooks for different steps
-   Low Priority
    -   GRN network inference
