# Background jobs
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(BSgenome.Nfurzeri.NCBI.Nfu20140520.custom)
  library(tidyverse)
  library(dtwclust)
})
source("C:/Users/Christopher He/OneDrive/UCSF/Singh Lab/diapause-sc-multiome/code/helper_functions.R")
setwd("C:/Users/Christopher He/OneDrive/UCSF/Singh Lab/diapause-sc-multiome/code/")

# so_neuron <- readRDS("../results_all/so_neurons.rds")
# so_neuron <- LinkPeaks(
#   object = so_neuron,
#   peak.assay = "ATAC",
#   expression.assay = "RNA"
# )
# saveRDS(so_neuron, "../results_all/so_neurons_links.rds")
# linkages <- as.data.frame(Links(so_neuron))
# write.csv(linkages, "../results_scratch/neuron_linkages.csv", quote = F)
################################################################
so_somite <- readRDS("../results_all/lineage_subsets/so_somite.rds")
so_somite <- LinkPeaks(
  object = so_somite,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)
saveRDS(so_somite, "../results_all/so_somite.rds")
linkages <- as.data.frame(Links(so_somite))
write.csv(linkages, "../results_scratch/linkages_somite.csv", quote = F)
################################################################
so_epidermal <- readRDS("../results_all/lineage_subsets/so_epidermal.rds")
so_epidermal <- LinkPeaks(
  object = so_epidermal,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)
saveRDS(so_epidermal, "../results_all/so_epidermal.rds")
linkages <- as.data.frame(Links(so_epidermal))
write.csv(linkages, "../results_scratch/linkages_epidermal.csv", quote = F)
################################################################
so_hatch <- readRDS("../results_all/lineage_subsets/so_hatch.rds")
so_hatch <- LinkPeaks(
  object = so_hatch,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)
saveRDS(so_hatch, "../results_all/so_hatch.rds")
linkages <- as.data.frame(Links(so_hatch))
write.csv(linkages, "../results_scratch/linkages_hatch.csv", quote = F)
################################################################
so_endothelial <- readRDS("../results_all/lineage_subsets/so_endothelial.rds")
so_endothelial <- LinkPeaks(
  object = so_endothelial,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)
saveRDS(so_endothelial, "../results_all/so_endothelial.rds")
linkages <- as.data.frame(Links(so_endothelial))
write.csv(linkages, "../results_scratch/linkages_endothelial.csv", quote = F)
################################################################
so_blood <- readRDS("../results_all/lineage_subsets/so_blood.rds")
so_blood <- LinkPeaks(
  object = so_blood,
  peak.assay = "ATAC",
  expression.assay = "RNA"
)
saveRDS(so_blood, "../results_all/so_blood.rds")
linkages <- as.data.frame(Links(so_blood))
write.csv(linkages, "../results_scratch/linkages_blood.csv", quote = F)