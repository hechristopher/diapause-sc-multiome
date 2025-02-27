# Background jobs
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(BSgenome.Nfurzeri.NCBI.Nfu20140520.custom)
  library(tidyverse) 
})


# build a joint neighbor graph using both assays
so <- FindMultiModalNeighbors(
  object = so,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 1:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
so <- RunUMAP(
  object = so,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "umap_wnn",
  verbose = TRUE
)

DimPlot(so, reduction = "umap_wnn") + NoLegend()

saveRDS(so, "../results_all/so_wnn.rds")