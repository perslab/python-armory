# Title: Functions to interact, plot and excel export data from Seurat objects.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017


####################################################################################################
########################################### WORKFLOW ###############################################
####################################################################################################

### Inside your script:
# xlsx.workbook <- createWorkbook(creator="PTimshel") # start excel
# .... call functions that adds sheets to your xlsx.workbook, e.g de_analysis()
# saveWorkbook(wb, file = "createWorkbookExample.xlsx", overwrite = TRUE) # write file



####################################################################################################
############################################# SETUP ################################################
####################################################################################################

library(openxlsx) # better than xlsx library because it does not depend on Java (which causes "java.lang.OutOfMemoryError" errors). Also openxlsx is being actively developed

library(dplyr)

library(gProfileR)


####################################################################################################
############################################# UTILS ################################################
####################################################################################################

remove_existing_excel_sheet <- function(excel_wb, sheet_name) {
  if (any(sheet_name %in% names(excel_wb))) { # check for existing sheet name
    # names(wb): get or set worksheet names for a Workbook object.
    # To avoid ERROR: Sheet '<sheet_name>' does not exist.
    print(sprintf("Sheet name %s already exists in the Excel workbook. It will be replaced", sheet_name))
    removeWorksheet(wb=excel_wb, sheet=sheet_name) # sheet should be a name or index of a worksheet
  }
}



add_cell_type_annotation <- function(df.in, df.cluster_annotation, cluster="cluster") {
  ### INPUT
  # df.in                   input data frame where column will be added 
  # df.cluster_annotation:  data frame with columns "cluster" and "annotation". 
  #                         Each row should map the cluster number ("cluster") to a cell type ("annotation").
  #                         *Only 1-1 mappings are allowed*
  # cluster                 string that specifies the column name of the cluster id column in df.in
  ### OUTPUT
  # a data frame with the column "annotation" added
  
  df.cluster_annotation <- as.data.frame(df.cluster_annotation) # in case we get a tibble, we convert it to a data frame, so we can use the df[, "x"] notation.
  df.in <- as.data.frame(df.in) # in case we get a tibble, we convert it to a data frame, so we can use the df[, "x"] notation.
  
  # Check that "cluster" and "annotation" values are unique (ensuring 1-1 mapping)
  if (anyDuplicated(df.cluster_annotation[,"cluster"]) |  anyDuplicated(df.cluster_annotation[,"annotation"]) ) { # anyDuplicated returns zero if there is NO duplicates
    stop("Error in add_cell_type_annotation: df.cluster_annotation contains non unique mappings. Please fix this.")
  }
  
  df.in$annotation <- as.character(plyr::mapvalues(x=df.in[,cluster], from=df.cluster_annotation[,"cluster"], to=as.character(df.cluster_annotation[,"annotation"]))) # add cell type information | as.character() avoids factor creation...
   # ^ notice that we do not need to load the whole plyr package - we can just refer directly to it.
  
  return(df.in)
}


####################################################################################################
######################################### PURE toExcel #############################################
####################################################################################################
# These functions primary purpose is to export data to Excel and does not return anything really useful.


################### DE GENES / CLUSTER MARKERS ##############################

df_to_excel <- function(df, excel_wb, sheet_name) {
  # Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addWorksheet(wb=excel_wb, sheetName=sheet_name)
  writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df), colNames=T, rowNames=F)
}

################### DE GENES / CLUSTER MARKERS ##############################

de_genes <- function(df.de_genes, df.cluster_annotation=NULL, colname_cluster="cluster",
                     excel_wb, sheet_name="de_genesX") {
  ### SEE de_genes_wide for documentation
  
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.de_genes <- add_cell_type_annotation(df.de_genes, df.cluster_annotation, cluster=colname_cluster)
  }
  
  # Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addWorksheet(wb=excel_wb, sheetName=sheet_name)
  writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.de_genes), colNames=T, rowNames=F)
  
}


################### TOP DE GENES / CLUSTER MARKERS: wide format ##############################

de_genes_wide <- function(df.de_genes, n_top_genes, df.cluster_annotation=NULL, colname_cluster="cluster",
                                 excel_wb, sheet_name="de_genesX.top") {
  ### INPUT
  # df.de_genes:              output data frame from Seurat FindAllMarkers() or FindMarkers().
  #                           data frame MUST contain columns "p_val" and "avg_logFC". Preferably also contain a column gene with genes (and not as rownames)
  #                           you may pre-filter the data frame to p_val_adj < 0.05.
  # n_top_genes:              number of top markers to write per cluster. 
  #                           Valid arguments: integer or 'max'. 
  #                           'max' will export the maximal number of markers, which is the number of marker genes of the cluster with the least marker genes.
  # cluster_annotation:       data frame with mapping from Seurat cluster id to cell type annotations (UNIQUE).
  #                           if this argument is provided, column name suffix will be the annotations and not cluster id.
  
  
  # SNIPPET df.de_genes (output from FindAllMarkers(), includes "gene" column
  # p_val  avg_logFC pct.1 pct.2    p_val_adj cluster    gene
  # ODC1    0.000000e+00 -0.8609366 0.069 0.275 0.000000e+00       0    ODC1
  # LTB     4.433493e-13  2.5034895 1.000 0.333 1.019703e-10       0     LTB
  # LDHB    1.153312e-08  1.8824861 0.552 0.451 2.652618e-06       0    LDHB
  
  group_cluster_var_name = "cluster" # default name for cluster id column from Seurat FindAllMarkers()
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.de_genes <- add_cell_type_annotation(df.de_genes, df.cluster_annotation, cluster=colname_cluster)
    group_cluster_var_name = "annotation" # add_cell_type_annotation() adds annotation column
  }
  
  n_max_genes <- df.de_genes %>% count(cluster) %>% pull(n) %>% min() # we need to find the least number of marker genes
  if (n_top_genes == "max") {
    n_top_genes <- n_max_genes
    print(sprintf("Max number of top markers selected: %s", n_top_genes))
  } else if (n_top_genes > n_max_genes) {
    stop(sprintf("n_top_genes is greater than the maximal allowed number of markers of %s", n_max_genes))
  }
  
  ### Top markers
  df.de_genes.top <- df.de_genes %>%
    group_by(UQ(rlang::sym(group_cluster_var_name))) %>% # group
    arrange(p_val, desc(avg_logFC)) %>% # sort
    slice(1:n_top_genes) # extract top markers
  
  
  list.tmp_res <- list()
  for (i in sort(unique(df.de_genes.top %>% pull(UQ(rlang::sym(group_cluster_var_name)))))) {
    df.tmp_cluster <- df.de_genes.top %>% filter(UQ(rlang::sym(group_cluster_var_name))==i)
    # df.tmp_cluster <- df.tmp_cluster %>% rename(pct_expr=pct.1, pct_expr_all_other_clusters=pct.2) # rename
    colnames(df.tmp_cluster) <- paste(colnames(df.tmp_cluster), sprintf("cluster.%s", i), sep=".")
    list.tmp_res[[i]] <- df.tmp_cluster
  }
  df.de_genes.top.wide <- bind_cols(list.tmp_res)
  
  # Write to excel file
  remove_existing_excel_sheet(excel_wb, sheet_name)
  addWorksheet(wb=excel_wb, sheetName=sheet_name)
  writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.de_genes.top.wide), colNames=T, rowNames=F)
  
}




###################################################################################################
####################################### COMPUTE FUNCTIONS #########################################
####################################################################################################
# These functions mostly compute something. Either as a wrapper around Seurat functions or as a new functionallity
# They also allow for expor to Excel



###################### Average Expression #####################

average_expression <- function(seurat_obj, colname_cluster_ident, log.space=T,
                               do.excel_export=F, excel_wb=NULL, sheet_name="avg_expr") {
  ### INPUT
  # colname_cluster_ident:                name of column of seurat @meta.data to group by, e.g. 'res.0.8' or 'annotation' or 'cell_types'.
  #                                       for each variable value in colname_cluster_ident and expression average will be calculated.
  # log.space                             if true, then return average in log-space (but averaging is done in non-log space). 
  #                                       returning average in log-space was the default behaviour of AverageExpression prior to Seurat 2.2.
  ### Output
  # data frame                            Averaged expression data. 

  df.avg_expr <- AverageExpression(SetAllIdent(seurat_obj, id=colname_cluster_ident), genes.use=NULL, return.seurat=FALSE, add.ident=NULL, use.scale=FALSE, use.raw=FALSE)
  # AverageDetectionRate()/AverageExpression() # Returns a data frame with genes as rows, identity classes as columns.
  if (log.space) {
    df.avg_expr <- log1p(df.avg_expr)
    sheet_name <- sprintf("%s.logspace", sheet_name) # update sheet name
  }
  df.avg_expr <- df.avg_expr %>% rownames_to_column(var="gene") # adding rownames as column.
  
  # NEW IN SEURAT VERSION 2.2: Output is in log-space when return.seurat = TRUE, otherwise it's in non-log space. Averaging is done in non-log space.
  # Use this to return average in log-space: log1p(AverageExpression(t.cells, show.progress = FALSE))
  # PRE-Version 2.2: "Output is in log-space, but averaging is done in non-log space."
  
  ### ABOUT the @data slot: normalized and log-transformed single cell expression but *NOT regressed*
  # We use the data in @data slot to calculate the average expression.
  # The data slot (object@data) stores normalized and log-transformed single cell expression. This maintains the relative abundance levels of all genes, and contains only zeros or positive values. See ?NormalizeData for more information.
  # This data is used for visualizations, such as violin and feature plots, most differential expression tests, finding high-variance genes, and as input to ScaleData (see below).
  
  
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.avg_expr), colNames=T, rowNames=F)
  }
  
  return(df.avg_expr)
}

###################### Average Dectection #####################

average_dectection_rate <- function(seurat_obj, colname_cluster_ident, 
                                    do.excel_export=F, excel_wb=NULL, sheet_name="avg_detection_rate") {
  # SEE average_expression
  
  # AverageDetectionRate()/AverageExpression() # Returns a matrix with genes as rows, identity classes as columns.
  df.avg_dectection_rate <- AverageDetectionRate(SetAllIdent(seurat_obj, id=colname_cluster_ident)) %>% rownames_to_column(var="gene")
  
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.avg_dectection_rate), colNames=T, rowNames=F)
  }
  
  return(df.avg_dectection_rate)
}



###################### DE ANALYSIS: ident specific (samples, groups, subpopulations) #####################


de_analysis <- function(seurat_obj, colname_ident, ident.1, ident.2, 
                                statistical_tests=c("wilcox", "negbinom"), 
                                do.excel_export=F, excel_wb=NULL, sheet_name=NULL) {
  ### INPUT
  # colname_ident:      name of any column in @meta.data that defines the grouping of ident.1 and ident.2.
  # ident.1 / ident.2   a string or vector specifying the grouping. E.g. ident.1="beta", ident.1=c(1,4) or ident.1=c("beta", "alpha")
  # test.use            a string or character vector specifying the statistical tests to run.
  ### OUTPUT
  # df.de_analysis:     data frame result from FindMarkers
  # adds the xlsx sheet to the specified workbook.
  

  if (is.null(sheet_name)) { # create sheet_name string if it is not provided
    sheet_name <- sprintf("de.%s.%s-vs-%s", colname_ident, paste0(ident.1, collapse="_"), paste0(ident.2, collapse="_"))
      # ^ paste0 makes the code work when ident.1 / ident.2 are vectors with multiple values
    print(sprintf("Writing to sheet name: %s", sheet_name))
  }
  
  seurat_obj <- SetAllIdent(seurat_obj, id=colname_ident) # set ident so we get the correct grouping
  
  list.findmarkers <- list()
  for (test.use in statistical_tests) {
    print(sprintf("Running FindMarkers for statistical test %s", test.use))
    df.res <- FindMarkers(seurat_obj, ident.1=ident.1, ident.2=ident.2, test.use=test.use)
    ### FindMarkers() OUT SNIPPET [data frame]. NOTE THAT GENE NAMERS ARE ROWNAMES AND NOT A COLUMN.
    # p_val avg_logFC pct.1 pct.2    p_val_adj
    # GNG11     8.794299e-19  6.014094   1.0 0.000 2.022689e-16
    # CLU       8.794299e-19  5.869606   1.0 0.000 2.022689e-16
    list.findmarkers[[test.use]] <- df.res %>% 
      rename_all(funs(paste0(names(df.res), ".", test.use))) %>% 
      rownames_to_column(var="gene") # add suffix to all columns and add gene as column
  }
  df.de_analysis <- plyr::join_all(list.findmarkers, type="full", by="gene") # plyr join list of data frame.
  

  df.de_analysis <- df.de_analysis %>% 
    rename( UQ(rlang::sym("logFC")) := UQ(rlang::sym(sprintf("avg_logFC.%s",statistical_tests[1]))) ) %>%  # renaming one of the avg_logFC columns. They are all the same, since they are not dependent on the statistical test.
    rename( UQ(rlang::sym("pct1")) := UQ(rlang::sym(sprintf("pct.1.%s",statistical_tests[1]))) ) %>% # renaming one pct.1 columns
    rename( UQ(rlang::sym("pct2")) := UQ(rlang::sym(sprintf("pct.2.%s",statistical_tests[1]))) ) %>% # renaming one pct.2 columns
    select(-starts_with("pct."), -starts_with("avg_logFC")) %>% # remove redundant columns
    arrange( UQ(rlang::sym(sprintf("p_val.%s", tail(statistical_tests, n=1)))) ) # sort by p_val for the last run of statistical_tests

  ## alternative renaming - WORKS, but clunky...
  # rename( UQ(rlang::sym(sprintf("logFC_%s_vs_%s", ident.1, ident.2))) := UQ(rlang::sym(sprintf("avg_logFC.%s",statistical_tests[1]))) ) %>%  # renaming one of the avg_logFC columns. They are all the same, since they are not dependent on the statistical test.
  # rename( UQ(rlang::sym(sprintf("pct_%s", ident.1))) := UQ(rlang::sym(sprintf("pct.1.%s",statistical_tests[1]))) ) %>% # renaming one pct.1 columns
  # rename( UQ(rlang::sym(sprintf("pct_%s", ident.2))) := UQ(rlang::sym(sprintf("pct.2.%s",statistical_tests[1]))) ) %>% # renaming one pct.2 columns
  
    
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.de_analysis), colNames=T, rowNames=F)
  }
  
  return(df.de_analysis)
}



############################## PER CLUSTER SAMPLE COMPOSITION ####################################

per_cluster_sample_composition <- function(seurat_obj, colname_cluster_ident, colname_group, df.cluster_annotation=NULL, colname_cluster="cluster",
                                           do.excel_export=F, excel_wb=NULL, sheet_name="clusters.composition") {
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
    # mutate(total_cells = sum(select(. , ends_with(".n")))) %>% # does not work
    arrange(as.numeric(cluster))
  # # A tibble: 19 x 5
  # # Groups:   cluster [19]
  # cluster hfd_ad_lib_lira.pct hfd_pair_fed.pct hfd_ad_lib_lira.n hfd_pair_fed.n
  # <chr>               <dbl>            <dbl>             <dbl>          <dbl>
  # 1       0                35.5             64.5               473            860
  # 2       1                35.0             65.0               457            850
  # 3       2                43.0             57.0               363            482
  
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.per_cluster_sample_composition.fmt <- add_cell_type_annotation(df.per_cluster_sample_composition.fmt, df.cluster_annotation, cluster=colname_cluster)
  }
  
  
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.per_cluster_sample_composition.fmt), colNames=T, rowNames=F)
  }
  
  return(df.per_cluster_sample_composition.fmt)
}



################### gProfileR enrichment ##############################

# ABOUT gProfileR Custom gene list as background (custom_bg): http://biit.cs.ut.ee/gprofiler/help.cgi?help_id=40
# In order to compute functional enrichments of gene lists, g:Profiler uses the backgound set of all organism-specific genes annotated in the Ensembl database.
# In several occasions, it is advisable to limit the background set for more accurate statistics. 
# For instance, one may use a custom background when the number of genes and corresponding probesets of a microarray platform is 
# considerably smaller than the number of known genes, or only genes of a specific chromosome are considered. 
# g:Profiler provides means to define the custom background as a mixed list of gene, probeset and protein IDs in the corresponding form field.
# It is also possible to select a predefined custom background from a list of popular microarray platforms. 

per_cluster_gprofiler_enrichment <- function(df.de_genes, organism, ordered_query=T, custom_bg="", 
                                         plot_fileout_prefix=NULL, 
                                         df.cluster_annotation=NULL, colname_cluster="cluster",
                                         do.excel_export=F, excel_wb=NULL, sheet_name=NULL) {
  ### SYNOPSIS
  # This function runs per cluster or cluster annotation.
  
  ### INPUT
  # df.de_genes:                          output data frame from Seurat differential expression tests (e.g. FindAllMarkers or FindMarkers)
  #                                 MUST contain columns "gene", "cluster" [or group_cluster_var_name], "avg_logFC", "p_val".
  #                                 should be FILTERED to your needs (e.g. p_val_adj < 0.05).
  # organism                        e.g. "hsapiens" or "mmusculus"
  # ordered_query                   for ranked gene lists to get GSEA style p-values.
  # custom_bg                       vector of gene names to use as a statistical background. 
  #                                 Should use the same gene identifiers (Ensembl IDs or gene names) as in the df.de_genes.
  #                                 default background is all annotated organism genes
  # plot_fileout_prefix             if argument specified, a png plot is generated with the provided argument value as prefix for the filename (full path, do NOT include file extension, e.g. png).
  #                                 the value of 'group_cluster_var_name' is used as suffix (e.g. cluster_1 or Oligodendrocytes_2)
  #                                 All directories will be created.
  #                                 EXAMPLE: /projects/timshel/sc-arc_lira/src/gprofiler_pngs/lira_vs_pf.cluster_markers.pval005
  
  ### OUTPUT
  # df.gprofiler                    data frame with results from gprofiler
  # .. to Excel                     writes results to Excel file.
  
  if (do.excel_export & is.null(sheet_name)) {
    stop("Please provide a sheet name or disable do.excel_export")
  }
  
  group_cluster_var_name = "cluster" # default name for cluster id column from Seurat FindAllMarkers()
  if (!is.null(df.cluster_annotation)) { # add cluster annotation if df.cluster_annotation is specified
    df.de_genes <- add_cell_type_annotation(df.de_genes, df.cluster_annotation, cluster=colname_cluster)
    group_cluster_var_name = "annotation" # add_cell_type_annotation() adds annotation column
  }
  
  # preparing data frame for gprofiler (sorting and grouping)
  df.de_genes <- df.de_genes %>%  
    group_by(UQ(rlang::sym(group_cluster_var_name))) %>%
    arrange(p_val, desc(avg_logFC)) # sort (this is needed if ordered_query=T)
    
  # run gprofiler - note that df.de_genes is ALREADY grouped, so we do not have to run group_by again.
  df.gprofiler <- df.de_genes %>% do(gprofiler(.$gene, organism=organism, ordered_query=ordered_query, significant=T, custom_bg=custom_bg))

  # make plot - if plot_fileout_prefix given
  # df.gprofiler_plots_return <- NULL # Default value
  if (!is.null(plot_fileout_prefix)) {
    dir.create(dirname(plot_fileout_prefix), recursive=T, showWarnings=F) # create outdir
    df.gprofiler_plots_return <- df.de_genes %>% do(data.frame(
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
  
  
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.gprofiler), colNames=T, rowNames=F)
  }
  
  return(df.gprofiler) 
  
}



gprofiler_enrichment <- function(df.de_genes, organism, ordered_query=T, custom_bg="",
                                 do.excel_export=F, excel_wb=NULL, sheet_name=NULL) {
  ### SYNOPSIS
  # This function runs a enrichment test on the genes provided in the "gene" column
  
  ### INPUT
  # SEE per_cluster_gprofiler_enrichment()
  
  # run gprofiler
  df.de_genes <- df.de_genes %>% 
    arrange(p_val, desc(avg_logFC)) # sort (this is needed if ordered_query=T)
    
  df.gprofiler <- gprofiler(df.de_genes$gene, organism=organism, ordered_query=ordered_query, significant=T, custom_bg=custom_bg)
  
  # Write to excel file
  if (do.excel_export) {
    remove_existing_excel_sheet(excel_wb, sheet_name)
    addWorksheet(wb=excel_wb, sheetName=sheet_name)
    writeData(wb=excel_wb, sheet=sheet_name, as.data.frame(df.gprofiler), colNames=T, rowNames=F)
  }
  
  return(df.gprofiler) 
}


####################################################################################################
################################### LEGACY CODE: xlsx library ######################################
####################################################################################################

### Inside your script:
# xlsx.workbook <- createWorkbook(type="xlsx") # start excel
# .... call functions that adds sheets to your xlsx.workbook, e.g de_analysis()
# saveWorkbook(xlsx.workbook, "1712XX_analysis.xlsx") # write file

### xlsx library
# ISSUE when writing a semi-large data frame to Excel: java.lang.OutOfMemoryError: GC overhead limit exceeded r xlsx
# REF: https://stackoverflow.com/a/27155704/6639640
# REF: https://stackoverflow.com/a/21937674/6639640
# As it is mentioned above run the options function at the beginning of your script before loading any libraries
# and if you are running it through Rstudio make sure you restart it before you run the script.
# options(java.parameters = "- Xmx8000m") # needed if writing semi-large data frame to Excel 
# library(xlsx)


# remove_existing_excel_sheet <- function(excel_wb, sheet_name) {
#   if (any(sheet_name %in% names(getSheets(excel_wb)))) { # check for existing sheet name
#     # getSheets returns a list of java object references each pointing to an worksheet. The list is named with the sheet names.
#     # To avoid ERROR: java.lang.IllegalArgumentException: The workbook already contains a sheet of this name
#     print(sprintf("Sheet name %s already exists in the Excel workbook. It will be replaced", sheet_name))
#     removeSheet(excel_wb, sheetName=sheet_name)
#   }
# }

# Write to excel file using xlsx library
# remove_existing_excel_sheet(excel_wb, sheet_name)
# addDataFrame(x=as.data.frame(df.de_genes.top.wide),
#              row.names=F, col.names=T,
#              sheet=createSheet(wb=excel_wb, sheetName=sheet_name))


