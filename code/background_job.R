# Background jobs
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(BSgenome.Nfurzeri.NCBI.Nfu20140520.custom)
  library(tidyverse) 
})
source("C:/Users/Christopher He/OneDrive/UCSF/Singh Lab/diapause-sc-multiome/code/helper_functions.R")

GlobalCelltypeDiaDE <- function(
    min_cell_number,
    timepoints_1 = c("dia6d", "dia1mo"),
    timepoints_2 = c("hb", "dev1d"),
    output_columns = c("pct.1", "pct.2","avg_log2FC", "p_val_adj"),
    downsample_to = F,
    assay = "RNA",
    seed = NULL
)
{
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  DefaultAssay(so) <- assay
  
  # get list of cell types which meet min cell number
  tp1_celltypes <- (as.data.frame(table((so@meta.data %>% filter(timepoint %in% timepoints_1))$annotation_broad)) %>%
                      filter(Freq > min_cell_number))$Var1
  tp2_celltypes <- (as.data.frame(table((so@meta.data %>% filter(timepoint %in% timepoints_2))$annotation_broad)) %>%
                      filter(Freq > min_cell_number))$Var1
  
  cell_types <- intersect(tp1_celltypes, tp2_celltypes)
  
  out <- data.frame(
    row.names = rownames(so),
    final_symbol = translate(rownames(so))
  )
  # global timepoint1 vs timepoint2
  cells.1 <- WhichCells(so, expression = timepoint %in% timepoints_1, seed = seed)
  cells.2 <- WhichCells(so, expression = timepoint %in% timepoints_2, seed = seed)
  
  if(downsample_to != F){
    cells.1 <- sample(cells.1, downsample_to)
    cells.2 <- sample(cells.2, downsample_to)
  }
  df <- FindMarkers(
    so,
    ident.1 = cells.1,
    ident.2 = cells.2,
    logfc.threshold = -Inf,
    min.pct = -Inf
  )
  df <- df[,output_columns]
  colnames(df) <- paste0(colnames(df), "_", "global")
  
  #merge with full out matrix
  out <- merge(out, df, by=0)
  rownames(out) <- out[,1]
  out <- out[,-1]
  
  
  #DE for each cell type
  for(i in 1:length(cell_types)){
    cell_type <- cell_types[i]
    
    #get barcodes for each timepoint (group)
    cells.1 = WhichCells(
      so,
      expression =
        (timepoint %in% timepoints_1) &
        (annotation_broad == cell_type),
      seed = seed
    )
    cells.2 = WhichCells(
      so,
      expression =
        (timepoint %in% timepoints_2) &
        (annotation_broad == cell_type),
      seed = seed
    )
    if(downsample_to != F){
      cells.1 <- sample(cells.1, downsample_to)
      cells.2 <- sample(cells.2, downsample_to)
    }
    
    #perform DE
    df <- FindMarkers(
      so,
      ident.1 = cells.1,
      ident.2 = cells.2,
      logfc.threshold = -Inf,
      min.pct = -Inf
    )
    df <- df[,output_columns]
    colnames(df) <- paste0(colnames(df), "_", cell_type)
    
    #merge with full out matrix
    out <- merge(out, df, by=0)
    rownames(out) <- out[,1]
    out <- out[,-1]
  }
  return(out)
}


so <- readRDS("../results_all/so_wnn.rds")


#iterative downsampling
#min_cells = 400, downsample to 300, iterate 100 times

n_iters <- 100
atac_iter <- vector('list', n_iters)
for(i in 1:n_iters){
  df <- GlobalCelltypeDiaDE(min_cell_number = 400, downsample_to = 300, assay = "ATAC")
  df$gene <- rownames(df)
  atac_iter[[i]] <- df
  print(paste0("iteration ", i, " complete"))
}

atac_iter <- bind_rows(atac_iter)

# Identify unique cell types by extracting column names
cell_types <- unique(gsub("(avg_log2FC|p_val_adj)_", "", grep("(avg_log2FC|p_val_adj)", names(atac_iter), value = TRUE)))

# Initialize an empty list to store summary dataframes
summary_list <- list()

# Iterate over each cell type and summarize data
for (cell_type in cell_types) {
  logFC_col <- paste0("avg_log2FC_", cell_type)
  pval_col <- paste0("p_val_adj_", cell_type)
  
  summary_df <- atac_iter %>%
    group_by(NCBI) %>%
    summarize(
      mean_logFC = mean(!!sym(logFC_col), na.rm = TRUE),
      freq_significant = mean(!!sym(pval_col) < 0.05, na.rm = TRUE)
    ) %>%
    rename_with(~ paste0(., "_", cell_type), -NCBI)
  
  summary_list[[cell_type]] <- summary_df
  print(paste0(cell_type, " done"))
}

# Combine all summaries into a single dataframe
final_atac_iter <- Reduce(function(x, y) full_join(x, y, by = "NCBI"), summary_list)

saveRDS(final_atac_iter, "../results_all/ct_DA_iter.rds")