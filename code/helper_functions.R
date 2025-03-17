# helper_functions.R

# translate gene orthologs between species
# acceptable column names: "N. furzeri (NCBI)", "N. furzeri Final Symbol", "
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
colors_broad <- c(
  "neuron_sema3c" = "#66FF00",
  "neuron_dscam" = "#C0FF00",
  "neuron" = "#32CD32",
  "differentiating_neuron" = "#ACE1AF",

  "neuron_progenitor" = "#13EEA4",
  "brain" = "#458B74", #ATAC ident
  "spinal_cord" = "#6959CD",
  "floor_plate" = "#177245", #ATAC ident
  "telencephalon/optic_cup" = "#698B22",
  "midbrain_hindbrain_pax2/5_en1/2" = "#008B8B",
  "midbrain_hindbrain_wnt1/lmx1b" = "#8DEEEE",
  "midbrain_hindbrain" = "#008B8B", #ATAC ident
 
  "somite" = "#FF7256",
  "head_kidney" = "#BF5700",
  "myotome/muscle" = "#CD2626",
  "primitive_heart" = "#8B0000",
  "tail_bud" = "#FFA54F",
  "neural_crest" = "#CDCD00",
  "pigment_cell" = "#CDCD00",

  "endothelial" = "#DAA520",
  "blood" = "#FFEC8B",
  "immune_cell" = "#8B814C",

  "epidermal" = "#F72585",
  "periderm" = "#7209B7",
  "otic_vesicle" = "#3A0CA3",
  "primitive_gut" = "#FF91AF",
  "pronephros" = "#C21E56",

  "fin_bud" = "#FFDEAD",
  "notochord" = "#355CBC",
  "hatching_gland" = "#536878",
  "hibernating_unknown" = "#C0C0C0"
)


# 
# names(colors_broad) <- idents
# colors_broad["hatching_gland"] <- "slategrey"

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
