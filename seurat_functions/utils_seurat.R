# Title: Utilities related to Seurat objects.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017



library(dplyr)


####################################################################################################
########################################### UTILS ##################################################
####################################################################################################



filter_genes_in_seurat_data.dataframe <- function(df, seurat_obj, colname_gene, do.print=T) {
  ### AIM         filter data frame to contain only genes in seurat data.
  ### INPUT
  # ....
  ### OUTPUT
  # df            a data frame, filtered to contain only genes in the seurat data. All columns are returned.
  df <- as.data.frame(df) # make sure we start off with a data frame (and not a tibble).
  genes_seurat_data <- seurat_obj@data@Dimnames[[1]] # vector of genes
  bool.genes_in_seurat_data <- df[,colname_gene] %in% genes_seurat_data # boolean
  if (do.print) {
    genes_not_found <- df[!bool.genes_in_seurat_data, colname_gene]
    print(sprintf("Filtered data frame. Number of genes not found in Seurat data: %s", length(genes_not_found)))
    print(genes_not_found)
  }
  df.filtered <- df[bool.genes_in_seurat_data, ] # subset
  return(df.filtered)
}

filter_genes_in_seurat_data.vector <- function(x, seurat_obj, do.print=T) {
  # SEE filter_genes_in_seurat_data.dataframe.
  # This function takes a vector as input.
  genes_seurat_data <- seurat_obj@data@Dimnames[[1]] # vector of genes
  bool.genes_in_seurat_data <- x %in% genes_seurat_data # boolean
  if (do.print) {
    genes_not_found <- x[!bool.genes_in_seurat_data]
    print(sprintf("Filtered vector. Number of genes not found in Seurat data: %s", length(genes_not_found)))
    print(genes_not_found)
  }
  x.filtered <- x[bool.genes_in_seurat_data] # subset
  return(x.filtered)
}

