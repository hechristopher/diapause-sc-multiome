# diapause-sc-multiome

#### This project contains analysis of single cell multiomic data from multiple timepoints in killifish development.

| Sample   | Cellranger Cell Number | Filtered Cell Number | Barcode Suffix |
|----------|------------------------|----------------------|----------------|
| 25som_1  | 10,177                 |                      | -1             |
| hb_1     | 9,774                  |                      | -2             |
| hb_2     | 10,510                 |                      | -3             |
| dia1d_1  | 6,575                  |                      | -4             |
| dia2d_1  | 10,736                 |                      | -5             |
| dia6d_1  | 5,993                  |                      | -6             |
| dia6d_2  | 17,082                 |                      | -7             |
| dia1mo_1 | 12,245                 |                      | -8             |
| dia3mo_1 | 7,928                  |                      | -9             |
| dev1d_1  | 12,501                 |                      | -10            |
| dev1d_2  | 20,000                 |                      | -11            |
| dev1d_3  | 12,563                 |                      | -12            |

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
