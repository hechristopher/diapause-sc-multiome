# Background jobs

source("C:/Users/Christopher He/OneDrive/UCSF/Singh Lab/diapause-sc-multiome/code/helper_functions.R")

setwd("C:/Users/Christopher He/OneDrive/UCSF/Singh Lab/diapause-sc-multiome/code/")


#set output and figure directories
outs_dir <- "../results_all/lineage_subsets/"
fig_dir <- "../results_all/lineage_subsets/"

so <- readRDS("../results_all/so_wnn.rds")

subset_process_RNA_ATAC <- function(so, rna.cts){
  #subset cells based on rna annotation
  rna.cells <- WhichCells(so, expression = annotation_broad %in% rna.cts)
  so_subset <- subset(so, cells = rna.cells)
  
  #process RNA
  DefaultAssay(so_subset) <- "RNA"
  so_subset <- NormalizeData(so_subset)
  so_subset <- FindVariableFeatures(so_subset)
  so_subset <- ScaleData(so_subset)
  so_subset <- RunPCA(so_subset, verbose = FALSE)
  so_subset <- FindNeighbors(object = so_subset, dims = 1:30)
  so_subset <- FindClusters(object = so_subset, resolution = 0.5)
  so_subset <- RunUMAP(so_subset, dims = 1:30, reduction.name = "umap.rna")
  
  #process ATAC
  DefaultAssay(so_subset) <- "ATAC"
  so_subset <- RunTFIDF(so_subset)
  so_subset <- FindTopFeatures(so_subset, min.cutoff = 'q0')
  so_subset <- RunSVD(so_subset)
  so_subset <- FindNeighbors(object = so_subset, dims = 1:30)
  so_subset <- FindClusters(object = so_subset, resolution = 0.5)
  so_subset <- RunUMAP(so_subset, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac") # optional
  
  #WNN
  so_subset <- FindMultiModalNeighbors(so_subset, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
  so_subset <- RunUMAP(so_subset, nn.name = "weighted.nn", reduction.name = "wnn.umap") # optional
  
  return(so_subset)
}

subset_list <- list(
  "1_neuron" = list("neuron", "neuron_progenitor", "differentiating_neuron"),
  "2_midbrain_hindbrain_pax_en" = list("midbrain_hindbrain_pax2/5_en1/2"),
  "3_telencephalon" = list("telencephalon/optic_cup"),
  "4_epidermal" = list("epidermal"),
  "5_periderm" = list("periderm"),
  "6_myotome" = list("myotome/muscle"),
  "7_somite" = list("somite"),
  "8_endothelial" = list("endothelial"),
  "9_blood" = list("blood"),
  "10_hatching_gland" = list("hatching_gland")
)


for(i in 1:length(subset_list)){
  name = names(subset_list)[i]
  rna.cts = subset_list[[i]]
  so_subset <- subset_process_RNA_ATAC(so, rna.cts)
  p1 <- DimPlot(so_subset, reduction = "umap.rna", group.by = "timepoint_group", cols = colors_timepoints)
  p2 <- DimPlot(so_subset, reduction = "umap.rna", group.by = "annotation_broad", cols = colors_broad, label = T)
  p3 <- DimPlot(so_subset, reduction = "umap.rna", group.by = "RNA_snn_res.0.5", label = T)
  
  
  p4 <- DimPlot(so_subset, reduction = "umap.atac", group.by = "timepoint_group", cols = colors_timepoints)
  p5 <- DimPlot(so_subset, reduction = "umap.atac", group.by = "annotation_broad", cols = colors_broad, label = T)
  p6 <- DimPlot(so_subset, reduction = "umap.atac", group.by = "ATAC_snn_res.0.5", label = T)
  
  p7 <- DimPlot(so_subset, reduction = "wnn.umap", group.by = "timepoint_group", cols = colors_timepoints)
  p8 <- DimPlot(so_subset, reduction = "wnn.umap", group.by = "annotation_broad", cols = colors_broad, label = T)
  
  SaveFigure((p1|p2|p3)/(p4|p5|p6)/(p7|p8) & NoLegend() & NoAxes(), paste0(name, "_subset_umaps"), width = 10, height = 10)
  saveRDS(so_subset, paste0(outs_dir, name, "_so.rds"))
}
