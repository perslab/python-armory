# Title: Functions to interact, plot and excel export data from Seurat objects.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017



library(dplyr)
library(ggplot2)





####################################################################################################
############################## GENE LIST FeaturePlot (per gene) ####################################
####################################################################################################

plot.gene_list_featureplot_per_gene <- function(seurat_obj, df.gene_list, col_name_lists="cell_type", col_name_genes="gene_name", dir_main_out="plots") {
  ### INPUT
  # df.gene_list:     data frame with two columns specified by col_name_lists and col_name_genes.
  ### OUTPUT
    # <dir_main_out>/<gene_list_name>/feature_plot.<gene_name>.pdf , e.g. plots.gene_lists/beta_cells/feature_plot.INS.pdf
  
  for (i in 1:nrow(df.gene_list)) {
    ### Constant
    # outdir <- "plots.gene_lists.feature" # outputs in currrent dir
    
    df.gene_list <- as.data.frame(df.gene_list) # # as.data.frame() to avoid error if a tibble is passed. 
      # tibbles does not allow for indexing via df[i,j]
      # ERROR : invalid or not-yet-implemented 'Matrix' subsetting)
    
    gene_list_name <- make.names(df.gene_list[i, col_name_lists]) # could also drop make.names()?
    gene_name <- df.gene_list[i, col_name_genes]
    tmp.outdir <- sprintf("%s/%s", dir_main_out, gene_list_name) # e.g. plots.gene_lists/beta_cells
    tmp.outfile <- sprintf("%s/feature_plot.%s.pdf", tmp.outdir, gene_name)
    
    if (!gene_name %in% rownames(seurat_obj@data)) { # check that gene is in the 10x/Seurat genes
      next
      print(sprintf("Gene name %s does not exist in Seurat data. Skipping it", gene_name))
    } else if (file.exists(tmp.outfile)) { # check if file exists
      next
    }
    
    tryCatch({
      print(sprintf("Processing #%s/#%s: %s | gene=%s", i, nrow(df.gene_list), tmp.outdir, gene_name))
      dir.create(tmp.outdir, recursive=TRUE, showWarnings=FALSE) # dir.create() does NOT crash if the directory already exists, it just prints out a warning.
      
      tmp.plot <- FeaturePlot(seurat_obj, gene_name, cols.use = c("lightgrey","blue"), do.return=T,  no.legend=F)
      ggsave(filename=tmp.outfile, w=8, h=6)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  } 
}

####################################################################################################
################## Dotplot AND FeaturePlot - TOP *CLUSTER* MARKER GENES ############################
####################################################################################################


plot.cluster_marker_genes <- function(seurat_obj, df.cluster_markers, col_cluster_id, n_top_markers=25, dir_main_out="plots") {
  ### INPUT
  # df.cluster_markers:             output data frame from Seurat FindAllMarkers(). MUST contain columns "cluster", "avg_logFC", "p_val"
  #                                 you should NOT pre-filter the data frame - it will be filtered in this function
  # col_cluster_id:     name of column in @meta.data to use for cluster id, e.g. "res.0.8"
  ### OUTPUT
  # <dir_main_out>/cluster_markers/feature_plot.<gene_name>.pdf , e.g. plots.gene_lists/beta_cells/feature_plot.INS.pdf
  
  loop.clusters <- unique(sort(as.numeric(seurat_obj@meta.data[,col_cluster_id])))
  for (i in loop.clusters) {
    ### Constant
    outdir <- "cluster_markers"
    
    tmp.outdir <- sprintf("%s/%s", dir_main_out, outdir)
    tmp.outfile.dotplot <- sprintf("%s/dotplot.cluster_%s.pdf", tmp.outdir, i)
    tmp.outfile.feature_plot <- sprintf("%s/feature_plot.cluster_%s.pdf", tmp.outdir, i)
    
    genes <- df.cluster_markers %>%
      filter(cluster==i, avg_logFC > 0) %>% # selecting cluster; positive marker
      arrange(p_val, desc(avg_logFC)) %>% # sorting
      slice(1:n_top_markers) %>% # select top genes
      pull(gene) # returns a vector of gene names
    
    print(sprintf("Processing #%s/#%s | cluster=%s", i, length(unique(loop.clusters)), i))
    dir.create(tmp.outdir, recursive=TRUE, showWarnings=FALSE) # dir.create() does NOT crash if the directory already exists, it just prints out a warning.
    
    ### DotPlot
    tmp.plot <- DotPlot(seurat_obj, genes, do.return=T, plot.legend=T, x.lab.rot=T)
    tmp.plot <- tmp.plot + labs(title=sprintf("cluster %s", i))
    ggsave(filename=tmp.outfile.dotplot, w=16, h=8)
    # --> w=20, h=12 is a good size if doing 20 plots (4 cols, 5 rows)
    
    ### FeaturPlot
    tmp.plot <- FeaturePlot(seurat_obj, genes, cols.use = c("lightgrey","blue"), do.return=T,  no.legend=F)
    ggsave(filename=tmp.outfile.feature_plot, w=20, h=12)
    # --> w=20, h=12 is a good size if doing 20 plots (4 cols, 5 rows)
    
  }
}  
  
