# helper_functions.R
suppressPackageStartupMessages({
  library(GOstats)
  library(GSEABase)
  library(Seurat)
  library(Signac)
  library(BSgenome.Nfurzeri.NCBI.Nfu20140520.custom)
  library(tidyverse)
})


# translate gene orthologs between species
# acceptable column names: "N. furzeri (NCBI)", "N. furzeri Final Symbol", "
translate <- function(genes, from = "N. furzeri (NCBI)", to = "N. furzeri Final Symbol", 
                      multiple = c("first", "all"), unlist_result = FALSE) {
  multiple <- match.arg(multiple)  # validate 'multiple' argument
  
  orthos <- readxl::read_excel("C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\gene_ortholog_translation_data.xlsx")
  
  # Remove NA from 'from' column
  orthos <- orthos[!is.na(orthos[[from]]), ]
  
  if (multiple == "first") {
    # Keep only first match per 'from'
    orthos_unique <- orthos[!duplicated(orthos[[from]]), ]
    row.names(orthos_unique) <- orthos_unique[[from]]
    translated <- orthos_unique[genes, ]
    result <- translated[[to]]
    return(result)
  }
  
  if (multiple == "all") {
    # Return all matches per input gene
    result <- lapply(genes, function(g) {
      matches <- orthos[orthos[[from]] == g, to]
      if (length(matches) == 0) NA else matches
    })
    names(result) <- genes
    
    if (unlist_result) {
      return(unlist(result, use.names = FALSE))  # flat unnamed vector
    } else {
      return(result)  # named list
    }
  }
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
                       "predia" = "#800080", "pre-dia" = "#800080",
                       "earlydia" = "#C6DBEF",
                       "middia" = "#6BAED6", "dia" = "#6BAED6",
                       "latedia" = "#08306B",
                       "dev" = "#FD8D3C"
                       )


#run go analysis on a set of genes

run_killifish_go <- function(genes, universe, ontolg = "BP", mingenes = 5, relenrich = 1, conditional = F, testDirection = "over", get_genes = F) {
  if(!file.exists("C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\diapause-sc-multiome\\GO\\gsc.rds")){
    frame = read.table(file ="C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\diapause-sc-multiome\\GO\\GO_killifish-human_best-hits_Nfur.txt",
                       header = T,
                       colClasses=c(rep("factor",3)))
    
    goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
    goFrame=AnnotationDbi::GOFrame(goframeData,organism="nfur")
    goAllFrame=AnnotationDbi::GOAllFrame(goFrame)
    gsc <- GSEABase::GeneSetCollection(goAllFrame, setType = GSEABase::GOCollection())
    saveRDS(gsc, "C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\diapause-sc-multiome\\GO\\gsc.rds")
  }
  gsc <- readRDS("C:\\Users\\Christopher He\\OneDrive\\UCSF\\Singh Lab\\diapause-sc-multiome\\GO\\gsc.rds")
  params <- Category::GSEAGOHyperGParams(
    name="My Custom GSEA based annot Params", 
    geneSetCollection=gsc, geneIds = genes, 
    universeGeneIds = universe, 
    ontology = ontolg,
    pvalueCutoff = 1,
    conditional = conditional, # Try True also here -- read documentation 
    testDirection = testDirection
  ) 
  # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over"
  print('starting HyperGTest')
  Over <- GOstats::hyperGTest(params)
  print("HyperGTest Done")
  enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))
  SummaryOver = data.frame(summary(Over), enrichment)
  FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]
  padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")
  FinalSummaryOver = data.frame(FilteredSummaryOver, padj)
  
  #GET GENES
  # isolate indexes for the go terms in final results
  ind.GO <- is.element(names(Over@goDag@nodeData@data), eval(parse(text=paste("FinalSummaryOver$", "GO",ontolg,"ID", sep=''))))
  selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]
  
  # get a go terms and genes in a new variable for all the terms in the results of enrichment
  goTerms <- lapply(selected.GO, 
                    function(x) x$geneIds)
  names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]
  
  # This will create a new file "genesForGOTerms.txt" that will have GO terms and genes in each gO terms
  # Genes can be duplicate, GO terms should not be
  # Number of Go terms or lines should be equal to the enriched go terms as in the other file generated by this script
  # This needs to be processed to generate the desired files
  
  if(nrow(FinalSummaryOver)==0){
    print("no GO Terms found")
    return(FinalSummaryOver)
  }
  genes_list <- list()
  for(i in 1:nrow(FinalSummaryOver)){
    goID <- FinalSummaryOver$GOBPID[i]
    goID_genes <- goTerms[[goID]]
    goID_genes <- goID_genes[goID_genes %in% genes]
    genes_list[[i]] <- goID_genes
  }
  FinalSummaryOver$genes <- genes_list
  
  return(FinalSummaryOver)
}
