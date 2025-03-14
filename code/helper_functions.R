# helper_functions.R

# translate gene orthologs between species
translate <- function(genes, from = "N. furzeri (NCBI)", to = "N. furzeri Final Symbol") {
  orthos <- readxl::read_excel("C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\gene_ortholog_translation_data.xlsx")
  
  orthos <- orthos[!duplicated(orthos[[from]]),]
  orthos <- orthos[!is.na(orthos[[from]]),]
  
  row.names(orthos) <- orthos[[from]]
  translated <- orthos[genes,]
  translated$genes <- genes
  translated <- translated %>% mutate({{ to }} := coalesce(.data[[to]], genes))
  return(translated[[to]])
}
# save figure
SaveFigure <- function(plots, path, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_dir, path, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_dir, path, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
  print(plots)
}

runSeuratRNAProcessing <- function(so, resolution = 0.8, ndims = 30){
  DefaultAssay(so) <- "RNA"
  so <- NormalizeData(so)
  so <- FindVariableFeatures(so)
  so <- ScaleData(so)
  so <- RunPCA(so)
  so <- RunUMAP(so, dims = 1:ndims)
  so <- FindNeighbors(so, dims = 1:ndims)
  so <- FindClusters(so, resolution = resolution)
  return(so)
}


colors_broad = Polychrome::palette36.colors(31)
colors_broad <- Polychrome::sortByHue(colors_broad)
idents <- c(
  "myotome/muscle",
  "head_kidney",
  "somite",
  "tail_bud",
  "endothelial",
  "blood",
  "neuron_sema3c",
  "neuron_dscam",
  "neuron",
  "differentiating_neuron",
  "spinal_cord",
  "floor_plate", #ATAC ident
  "neuron_progenitor",
  "brain", #ATAC ident
  "midbrain_hindbrain_pax2/5_en1/2",
  "midbrain_hindbrain_wnt1/lmx1b",
  "midbrain_hindbrain", #ATAC ident
  "telencephalon/optic_cup",
  "hibernating_unknown",
  "epidermal",
  "neural_crest",
  "periderm",
  "otic_vesicle",
  "primitive_gut",
  "pronephros",
  "hatching_gland",
  "pigment_cell",
  "notochord",
  "fin_bud",
  "primitive_heart",
  "immune_cell"
)
names(colors_broad) <- idents

colors_timepoints <- c("25som"= "#E496EF", "25som_1"= "#E496EF",
                       "hb"= "#800080", "hb_1"= "#800080", "hb_2"= "#800080",
                       "dia1d"= "#C6DBEF", "dia1d_1"= "#C6DBEF",
                       "dia2d"= "#C6DBEF", "dia2d_1"= "#C6DBEF",
                       "dia6d"= "#6BAED6", "dia6d_1"= "#6BAED6", "dia6d_2"= "#6BAED6",
                       "dia1mo"= "#2171B5",  "dia1mo_1"= "#2171B5",
                       "dia3mo"= "#08306B",  "dia3mo_1"= "#08306B",
                       "dev1d"= "#FD8D3C", "dev1d_1"= "#FD8D3C", "dev1d_2"= "#FD8D3C", "dev1d_3"= "#FD8D3C",
                       "predia" = "#800080",
                       "earlydia" = "#C6DBEF",
                       "middia" = "#6BAED6", 
                       "latedia" = "#08306B",
                       "dev" = "#FD8D3C"
                       )
