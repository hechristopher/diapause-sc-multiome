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

k_clusters <- 4
clusts0 <- tsclust(t(scaled_data), type='partitional',
                   k=2:10, # to check different values
                   distance='dtw_basic', centroid='pam', trace=T)