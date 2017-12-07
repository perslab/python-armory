# Title: Functions to interact, plot and excel export data from Seurat objects.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017


####################################################################################################
########################################### WORKFLOW ###############################################
####################################################################################################

### Inside your script:
# xlsx.workbook <- createWorkbook(type="xlsx") # start excel
# .... call functions that adds sheets to your xlsx.workbook, e.g toExcel.de_analysis()
# saveWorkbook(xlsx.workbook, "1712XX_analysis.xlsx") # write file


####################################################################################################
############################################# SETUP ################################################
####################################################################################################

library(xlsx)
library(dplyr)

library(gProfileR)


####################################################################################################
############################################# UTILS ################################################
####################################################################################################

remove_existing_excel_sheet <- function(excel_wb, sheet_name) {
  if (any(sheet_name %in% names(getSheets(excel_wb)))) { # check for existing sheet name
    # getSheets returns a list of java object references each pointing to an worksheet. The list is named with the sheet names.
    # To avoid ERROR: java.lang.IllegalArgumentException: The workbook already contains a sheet of this name
    print(sprintf("Sheet name %s already exists in the Excel workbook. It will be replaced", sheet_name))
    removeSheet(excel_wb, sheetName=sheet_name)
  }
}

add_cell_type_annotation <- function(df.in, df.cluster_annotation, colname_cluster_id="cluster") {
  ### INPUT
  # df.in                   input data frame where column will be added  
  # df.cluster_annotation:  data frame with columns "cluster_id" and "annotation"
  # colname_cluster_id      string that specifies the column name of the cluster id column in df.in
  ### OUTPUT
  # a data frame with the column "annotation" added
   df.in$annotation <- as.character(plyr::mapvalues(x=df.in[,colname_cluster_id], from=cluster_id, to=annotation)) # add cell type information | as.character() avoids factor creation...
   # ^ notice that we do not need to load the whole plyr package - we can just refer directly to it.
  return(df.in)
}


####################################################################################################
############################################ toExcel ###############################################
####################################################################################################

###################### DE ANALYSIS: ident specific (samples, groups, subpopulations) #####################


toExcel.de_analysis <- function(seurat_obj, colname_ident, ident.1, ident.2, 
                                test.use="wilcox", 
                                excel_wb, sheet_name=NULL) {
  ### INPUT
  # colname_ident:      name of any column in @meta.data that defines the grouping of ident.1 and ident.2.
  # ident.1 / ident.2   a string or vector specifying the grouping. E.g. ident.1="beta", ident.1=c(1,4) or ident.1=c("beta", "alpha")
  ### OUTPUT
  # df.de_analysis:     data frame result from FindMarkers
  # adds the xlsx sheet to the specified workbook.
  

  if (is.null(sheet_name)) { # create sheet_name string if it is not provided
    sheet_name <- sprintf("de.%s.%s-vs-%s", colname_ident, paste0(ident.1, collapse="_"), paste0(ident.2, collapse="_"))
      # ^ paste0 makes the code work when ident.1 / ident.2 are vectors with multiple values
    print(sprintf("Writing to sheet name: %s", sheet_name))
  }
  
  seurat_obj <- SetAllIdent(seurat_obj, id=colname_ident) # set ident so we get the correct grouping
  print(seurat_obj@ident %>% count())
  
  df.de_analysis <- FindMarkers(seurat_obj, ident.1=ident.1, ident.2=ident.2, test.use=test.use)
  
  ### Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addDataFrame(x=as.data.frame(df.de_analysis),
               row.names=F, col.names=T,
               sheet=createSheet(wb=excel_wb, sheetName=sheet_name))
  
  return(df.de_analysis)
}



############################## PER CLUSTER SAMPLE COMPOSITION ####################################

toExcel.per_cluster_sample_composition <- function(seurat_obj, colname_cluster_ident, colname_group, excel_wb, sheet_name="clusters.composition") {
  ### INPUT
  # colname_cluster_ident:                name of column of seurat @meta.data to group by, e.g. 'res.0.8' or 'annotation' or 'cell_types'
  # colname_group:                        grouping variable for the composition analysis, e.g 'sample' or 'group'
  # excel_wb:                              excel workbook object from the "xlsx" package
  ### OUTPUT
  # df.per_cluster_sample_composition.fmt   data frame
  # adds the xlsx sheet to the specified workbook (the modification to the workbook is done in place)
  
  ### Relative frequencies / proportions with dplyr
  # "When you group by multiple variables, each summary peels off one level of the grouping." 
  # --> the grouping order matters
  # REF: https://stackoverflow.com/a/24576703
  df.per_cluster_sample_composition <- seurat_obj@meta.data %>% 
    group_by(UQ(rlang::sym(colname_cluster_ident)), UQ(rlang::sym(colname_group))) %>%
    summarize(n=n()) %>%
    mutate(pct=round(n/sum(n)*100,1)) %>%
    rename(cluster=UQ(rlang::sym(colname_cluster_ident)))
  
  ### Making the format of the table nicer.
  # REF: https://stackoverflow.com/questions/30592094/r-spreading-multiple-columns-with-tidyr
  df.per_cluster_sample_composition.fmt <- df.per_cluster_sample_composition %>% 
    gather(tmp_key, tmp_value, n,pct) %>%
    unite(tmp_unite,UQ(rlang::sym(colname_group)),tmp_key, sep=".") %>%
    spread(tmp_unite, tmp_value, fill=0) %>% # fills 'missing combinations' with zero 
    select(-ends_with("n"), ends_with("n")) %>% # puts count ".n" clusters to the end
    arrange(as.numeric(cluster))


  ### Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addDataFrame(x=as.data.frame(df.per_cluster_sample_composition.fmt),
               row.names=F, col.names=T,
               sheet=createSheet(wb=excel_wb, sheetName=sheet_name))
  
  return(df.per_cluster_sample_composition.fmt)
}

################### CLUSTER MARKERS ##############################

toExcel.cluster_markers <- function(df.cluster_markers, excel_wb, sheet_name="clusters.markers", df.cluster_annotation=NULL) {
  ### SEE toExcel.cluster_markers_wide for documentation
  
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.cluster_markers <- add_cell_type_annotation(df.cluster_markers, df.cluster_annotation, colname_cluster_id="cluster")
  }
  
  ### Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addDataFrame(x=as.data.frame(df.cluster_markers),
               row.names=F, col.names=T,
               sheet=createSheet(wb=excel_wb, sheetName=sheet_name))
}


################### TOP *CLUSTER* MARKERS: long vs wide format ##############################

toExcel.cluster_markers_wide <- function(df.cluster_markers, n_top_markers, excel_wb, sheet_name="clusters.top_markers", df.cluster_annotation=NULL) {
  ### INPUT
  # df.cluster_markers:             output data frame from Seurat FindAllMarkers().
  #                                 you may pre-filter the data frame to p_val_adj < 0.05.
  # n_top_markers:                  number of top markers to write per cluster. 
  #                                 Valid arguments: integer or 'max'. 
  #                                   'max' will export the maximal number of markers, which is the number of marker genes of the cluster with the least marker genes.
  # cluster_annotation [optional]:  data frame with mapping from Seurat cluster id to cell type annotations (UNIQUE).
  #                                 if this argument is provided, column name suffix will be the annotations and not cluster id.
  
  group_cluster_var_name = "cluster" # default name for cluster id column from Seurat FindAllMarkers()
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.cluster_markers <- add_cell_type_annotation(df.cluster_markers, df.cluster_annotation, colname_cluster_id="cluster")
    group_cluster_var_name = "annotation" # add_cell_type_annotation() adds annotation column
  }
  
  n_max_markers <- df.cluster_markers %>% count(cluster) %>% pull(n) %>% min() # we need to find the least number of marker genes
  if (n_top_markers == "max") {
    n_top_markers <- n_max_markers
    print(sprintf("Max number of top markers selected: %s", n_top_markers))
  } else if (n_top_markers > n_max_markers) {
    stop(sprintf("n_top_markers is greater than the maximal allowed number of markers of %s", n_max_markers))
  }
  
  ### Top markers
  df.cluster_markers.top <- df.cluster_markers %>%
    group_by(UQ(rlang::sym(group_cluster_var_name))) %>% # group
    arrange(p_val, desc(avg_logFC)) %>% # sort
    slice(1:n_top_markers) # extract top markers
  
  
  list.tmp_res <- list()
  for (i in sort(unique(df.cluster_markers.top %>% pull(UQ(rlang::sym(group_cluster_var_name)))))) {
    df.tmp_cluster <- df.cluster_markers.top %>% filter(UQ(rlang::sym(group_cluster_var_name))==i)
    # df.tmp_cluster <- df.tmp_cluster %>% rename(pct_expr=pct.1, pct_expr_all_other_clusters=pct.2) # rename
    colnames(df.tmp_cluster) <- paste(colnames(df.tmp_cluster), sprintf("cluster.%s", i), sep=".")
    list.tmp_res[[i]] <- df.tmp_cluster
  }
  df.cluster_markers.top.wide <- bind_cols(list.tmp_res)
  
  ### Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addDataFrame(x=as.data.frame(df.cluster_markers.top.wide),
               row.names=F, col.names=T,
               sheet=createSheet(wb=excel_wb, sheetName=sheet_name))
}



################### gProfileR enrichment ##############################

# ABOUT Custom gene list as background (custom_bg): http://biit.cs.ut.ee/gprofiler/help.cgi?help_id=40
# In order to compute functional enrichments of gene lists, g:Profiler uses the backgound set of all organism-specific genes annotated in the Ensembl database.
# In several occasions, it is advisable to limit the background set for more accurate statistics. 
# For instance, one may use a custom background when the number of genes and corresponding probesets of a microarray platform is 
# considerably smaller than the number of known genes, or only genes of a specific chromosome are considered. 
# g:Profiler provides means to define the custom background as a mixed list of gene, probeset and protein IDs in the corresponding form field.
# It is also possible to select a predefined custom background from a list of popular microarray platforms. 

toExcel.gprofiler_enrichment <- function(df.de, organism, ordered_query=T, custom_bg="", 
                                         plot_fileout_prefix=NULL, 
                                         excel_wb, sheet_name, df.cluster_annotation=NULL) {
  ### INPUT
  # df.de:                          output data frame from Seurat differential expression tests (e.g. FindAllMarkers or FindMarkers)
  #                                 MUST contain columns "gene", "cluster" [or group_cluster_var_name], "avg_logFC", "p_val".
  #                                 should be FILTERED to your needs (e.g. p_val_adj < 0.05).
  # organism                        e.g. "hsapiens" or "mmusculus"
  # ordered_query                   for ranked gene lists to get GSEA style p-values.
  # custom_bg                       vector of gene names to use as a statistical background. 
  #                                 Should use the same gene identifiers (Ensembl IDs or gene names) as in the df.cluster_markers.
  #                                 default background is all annotated organism genes
  # plot_fileout_prefix             if argument specified, a png plot is generated with the provided argument value as prefix for the filename (full path, do NOT include file extension, e.g. png).
  #                                 the value of 'group_cluster_var_name' is used as suffix (e.g. cluster_1 or Oligodendrocytes_2)
  #                                 All directories will be created.
  #                                 EXAMPLE: /projects/timshel/sc-arc_lira/src/gprofiler_pngs/lira_vs_pf.cluster_markers.pval005
  
  ### OUTPUT
  # df.gprofiler                    data frame with results from gprofiler
  # .. to Excel                     writes results to Excel file.
  
  group_cluster_var_name = "cluster" # default name for cluster id column from Seurat FindAllMarkers()
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.cluster_markers <- add_cell_type_annotation(df.cluster_markers, df.cluster_annotation, colname_cluster_id="cluster")
    group_cluster_var_name = "annotation" # add_cell_type_annotation() adds annotation column
  }
  
  # preparing data frame for gprofiler (sorting and grouping)
  df.de <- df.de %>%  
    group_by(UQ(rlang::sym(group_cluster_var_name))) %>%
    arrange(p_val, desc(avg_logFC)) # sort (this is needed if ordered_query=T)
    
  # run gprofiler - note that df.de is ALREADY grouped, so we do not have to run group_by again.
  df.gprofiler <- df.de %>% do(gprofiler(.$gene, organism=organism, ordered_query=ordered_query, significant=T, custom_bg=custom_bg))

  # make plot - if plot_fileout_prefix given
  # df.gprofiler_plots_return <- NULL # Default value
  if (!is.null(plot_fileout_prefix)) {
    dir.create(dirname(plot_fileout_prefix), recursive=T, showWarnings=F) # create outdir
    df.gprofiler_plots_return <- df.de %>% do(data.frame(
      return_value_png_bool = gprofiler(.$gene, organism=organism, ordered_query=ordered_query, significant=T, custom_bg=custom_bg, 
                                                                          png_fn=sprintf("%s.%s.png", plot_fileout_prefix, unique(.$UQ(rlang::sym(group_cluster_var_name)))) ),
      png_filename = sprintf("%s.group_%s.png", plot_fileout_prefix, unique(.$UQ(rlang::sym(group_cluster_var_name))))
      #png_filename = sprintf("%s.%s.png", plot_fileout_prefix, unique(. %>% pull(UQ(rlang::sym(group_cluster_var_name))))) # DOES NOT WORK:  Error in unique.default(. %>% pull(cluster)) : unique() applies only to vectors 
      #png_filename1 = sprintf("%s.%s.png", plot_fileout_prefix, . %>% pull(UQ(rlang::sym(group_cluster_var_name))) %>% unique() ), # DOES NOT WORK:   Error in as.vector(x, "character") :  cannot coerce type 'closure' to vector of type 'character' 
      #png_filename2 = sprintf("%s.%s.png", plot_fileout_prefix, unique(.$cluster)), # WORKS - BUT not flexible to group_cluster_var_name
      ))
      # gprofiler() returns TRUE or FALSE when png_fn argument is given
      # EXPR in dplyr::do(EXPR) MUST return a data frame. Hence we wrap it in as.data.frame()
    print(df.gprofiler_plots_return) # a 3-column data frame (group_cluster_var_name, return_value_png_bool, png_filename)
  }
  
  
  ### Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addDataFrame(x=as.data.frame(df.gprofiler),
               row.names=F, col.names=T,
               sheet=createSheet(wb=excel_wb, sheetName=sheet_name))
  
  return(df.gprofiler) 
  
}








