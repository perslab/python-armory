# Title: Script to analyze and align two conditations directly from 'raw' 10x data.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017




######################################################################
############################## USAGE #################################
######################################################################

# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions lira_vs_pf |& tee analyze_two_conditions_from_raw.lira_vs_pf.out.txt
# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions hfd_vs_chow |& tee analyze_two_conditions_from_raw.hfd_vs_chow.out.txt
# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions hfd_vs_pf |& tee analyze_two_conditions_from_raw.hfd_vs_pf.out.txt

# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions hfd_vs_lira --n_cores 10 |& tee analyze_two_conditions_from_raw.hfd_vs_lira.out.txt
# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions pf_vs_chow --n_cores 10 |& tee analyze_two_conditions_from_raw.pf_vs_chow.out.txt


# time Rscript analyze_two_conditions_from_raw.R --analysis_conditions lira_vs_pf --output_prefix testX --test_run

######################################################################
########################### OptParse ################################
######################################################################
library(optparse)

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

option_list <- list( 
  make_option("--analysis_conditions", type="character",
              help = "Choose from 'lira_vs_pf', 'hfd_vs_chow', 'hfd_vs_pf', 'hfd_vs_lira', 'pf_vs_chow'"),
  make_option("--output_prefix", type="character", default=NULL,
              help = "Prefix for output files"),
  make_option("--test_run", action="store_true", default=FALSE, 
              help="Run in 'test mode': load only 500 cells and a limited number of highly expressed genes."),
  make_option("--n_cores", type="integer", default=10,
              help = "Number of cores to use for parallelization [default %default]")
  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


analysis_conditions <- opt$analysis_conditions # get argument


if (is.null(opt$output_prefix)) { # if no argument given, we set the output prefix to analysis_conditions
  output_prefix <- analysis_conditions
  print(sprintf("No output_prefix given. Using output_prefix = %s", output_prefix))
} else {
  output_prefix <- opt$output_prefix
  print(sprintf("Got output_prefix %s", output_prefix))
}


FLAG_TEST_RUN <- opt$test_run # set to TRUE to load only 500 cells and a limited number of highly expressed genes
# ^ uses Seurat.CreateObj.test_small_dataset()

N_CORES <- opt$n_cores


# IMPORTANT: condition1 and condition2 names *MUST* match the folder names in /data/sc-10x/data-runs/170612-perslab-arc_lira/
if (analysis_conditions == "lira_vs_pf") {
  condition1 <-"hfd_ad_lib_lira"
  condition2 <-"hfd_pair_fed"
} else if (analysis_conditions == "hfd_vs_chow") {
  condition1 <-"hfd_ad_lib"
  condition2 <-"chow"
} else if (analysis_conditions == "hfd_vs_pf") {
  condition1 <-"hfd_ad_lib"
  condition2 <-"hfd_pair_fed"
} else if (analysis_conditions == "hfd_vs_lira") {
  condition1 <-"hfd_ad_lib"
  condition2 <-"hfd_ad_lib_lira"
} else if (analysis_conditions == "pf_vs_chow") {
  condition1 <-"hfd_pair_fed"
  condition2 <-"chow"
} else {
  stop("Wrong argument for analysis_conditions")
}

######################################################################
########################### CONSTANTS ################################
######################################################################


print(sprintf("========== RUNNING analysis_conditions = %s ==========", analysis_conditions))


res.primary <- 0.8 # the primary resolution - will be used for cluster marker identification.
res2calculate <- c(0.4, 0.6, 0.8, 1.2, 1.6, 2, 4, 10) # resolutions to calculate for FindClusters
# OBS: res.primary must be included in res2calculate
# "name" in metadata is set via paste0("res.", r), so res=2.0 will be "res.2" and NOT "res.2.0"
# sprintf("res.%s", 2.0) --> "res.2"
# paste0("res.", 2.0) --> "res.2"


N_DIMS <- 20 # dims.use for CCA (and tSNE)
N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS <- 10 # used for condition_de_genes(seurat_obj)


######################################################################
######################### DENPENDENCIES ##############################
######################################################################

library(Seurat)
library(tidyverse)

library(parallel)

# tibble (part of tidyverse)
# plyr (plyr::join_all)


######################################################################
########################### FUNCTIONS ################################
######################################################################

Seurat.CreateObj.test_small_dataset <- function(list_name, list_of_matrix_data) {
  # ONLY USED FOR TESTING ON A SMALL DATASET
  
  ### METHOD 1: you know the number of genes and cells (but might be 'low quality' cells/genes)
  #N_GENES <- 2000; N_CELLS <- 500 # make sure to include enough genes for HVG/PCA/CCA/tSNE
  #data <- as.matrix(list_of_matrix_data[[list_name]])[1:N_GENES, 1:N_CELLS] # ROWS=GENES; COLUMNS=CELLS
  # ^ this way may not be the best way of subsetting if you the genes in the datasets are not in the same order.
  # ^ then you will end up with data sets with few overlapping genes, and the seurat test might not run optimally.
  #seurat_obj <- CreateSeuratObject(data, min.cells=0, min.genes=0, project=list_name) 
  # no filtering of cells or genes.
  # project specifies the 'orig.ident' column in @meta.data
  
  ### METHOD 2: you QC genes/cells and control the number of cells (SubsetData)
  seurat_obj <- CreateSeuratObject(raw.data=list_of_matrix_data[[list_name]], min.cells=300, min.genes=1500, project=list_name)
  # seurat_obj <- SubsetData(seurat_obj, cells.use=WhichCells())
  seurat_obj <- SubsetData(seurat_obj, max.cells.per.ident=500, random.seed=1) # alternative: SetAllIdent(seurat_obj, id="orig.ident")
  # ^ SubsetData calls WhichCells.
  # ^ Seurat sets 'project name' (orig.ident) as default identity
  # RESULT:
  # 6581 genes across 500 samples.
  # 8410 genes across 500 samples.
}


Seurat.CreateObj <- function(list_name, list_of_matrix_data) {
  seurat_obj <- CreateSeuratObject(raw.data=list_of_matrix_data[[list_name]], min.cells=3, min.genes=500, project=list_name) # project specifies the 'orig.ident' column in @meta.data
}


# Add mito/ribo meta data for each experiment
Seurat.AddMitoRiboPct <- function(seurat_obj) {
  ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat_obj@data), value = TRUE)
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_obj@data), value = TRUE, ignore.case=TRUE)
  percent.ribo <- Matrix::colSums(seurat_obj@raw.data[ribo.genes, ])/Matrix::colSums(seurat_obj@raw.data)
  percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.ribo, col.name = "percent.ribo")
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")
}

Seurat.Filter <- function(seurat_obj) {
  seurat_obj <- FilterCells(object = seurat_obj, subset.names = c("nGene", "percent.mito"), 
                            low.thresholds = c(200, -Inf), high.thresholds = c(4000,.3))
}

#Normalize,scale,find var genes
Seurat.NormAndScale <- function(seurat_obj) {
  seurat_obj<-NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
  seurat_obj<-ScaleData(seurat_obj, vars.to.regress = c("nUMI","percent.mito"))
  #seurat_obj<-ScaleData(seurat_obj, vars.to.regress = c("nUMI","percent.mito"), use.umi=T, model.use="negbinom") # too slow
  seurat_obj<-FindVariableGenes(seurat_obj, do.plot=F)
}


condition_de_genes <- function(seurat_obj) {
  # Find differentially eseurat_objpressed genes between the two conditions
  # *OBS* note that this function will return NULL if the 'if statement' is not true.
  
  statistical_tests <- c("wilcox", "negbinom")
  
  flag_condition1_ok <- sum(seurat_obj@meta.data$condition == condition1) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS
  flag_condition2_ok <- sum(seurat_obj@meta.data$condition == condition2) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS
  if (!(flag_condition1_ok & flag_condition2_ok)) { # if one of the flags are FALSE, the function will return NULL
    return(NULL)
  }
  
  
  list.findmarkers <- list()
  for (test.use in statistical_tests) {
    print(sprintf("Running FindMarkers for statistical test %s", test.use))
    df.res <- FindMarkers(seurat_obj, ident.1=condition1, ident.2=condition2, test.use=test.use) # Returns data frame
    # OBS: SetAllIdent(seurat_obj, "condition") MUST have been called on the Seurat object prior to this, or you will get an error "Identity : hfd_ad_lib_lira not found."
    
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
    rename( avg_logFC = logFC ) %>%  # renaming BACK AGAIN
    rename( UQ(rlang::sym("pct.1")) := pct1 ) %>%  # renaming BACK AGAIN, not sure if UQ(rlang::sym()) is needed
    rename( UQ(rlang::sym("pct.2")) := pct2 ) %>%  # renaming BACK AGAIN, not sure if UQ(rlang::sym()) is needed
    mutate(
      p_val = UQ(rlang::sym(sprintf("p_val.%s", tail(statistical_tests, n=1)))),
      p_val_adj = UQ(rlang::sym(sprintf("p_val_adj.%s", tail(statistical_tests, n=1))))
    ) %>%
  arrange( UQ(rlang::sym(sprintf("p_val.%s", tail(statistical_tests, n=1)))) ) # sort by p_val for the last run of statistical_tests
  
  ### SNITPPIT df.de_analysis
  # gene p_val.wilcox      logFC  pct1  pct2 p_val_adj.wilcox p_val.negbinom p_val_adj.negbinom
  # 1 FP671120.3 1.381887e-32 -1.2566769 0.979 1.000     4.588694e-28              0                  0
  # 2 FP236383.2 1.381887e-32 -1.2567393 0.979 1.000     4.588694e-28              0                  0
  # 3   FP236383 1.381887e-32 -1.2568003 0.979 1.000     4.588694e-28              0                  0
  return(df.de_analysis)
}




###########################################################################
############################# READ DATA ###################################
###########################################################################
print("Reading data")

matrix.dat1 <- Read10X(paste0("/raid5/data/sc-10x/data-runs/170612-perslab-arc_lira/",sprintf("agg-arc_lira-%s-no_norm", condition1),"/outs/filtered_gene_bc_matrices_mex/mm10/"))
matrix.dat2 <- Read10X(paste0("/raid5/data/sc-10x/data-runs/170612-perslab-arc_lira/",sprintf("agg-arc_lira-%s-no_norm", condition2),"/outs/filtered_gene_bc_matrices_mex/mm10/"))
# returns sparseMatrix. ROWS=GENES; COLUMNS=CELLS

list_of_matrix_data <-list(matrix.dat1, matrix.dat2) # *OBS*: the ORDER is important.
names(list_of_matrix_data) <- c(condition1, condition2) # IMPORTANT: *NAMED* LIST. Seurat.CreateObj() needs a named object
# names(list_of_matrix_data) ---> "hfd_ad_lib_lira" "hfd_pair_fed" 
# ^ the named list was not preserved in list_of_seurat_objs.

##############################################################################
############################# PROCESS DATA ###################################
##############################################################################
print("Processing data")


# Load data - create seurat objects for each experiment
if (FLAG_TEST_RUN) {
  print("***RUNNING TEST DATA MODE***")
  list_of_seurat_objs <- lapply(names(list_of_matrix_data), Seurat.CreateObj.test_small_dataset, list_of_matrix_data) # TEST on small data set
  print(list_of_seurat_objs)
} else {
  print("***RUNNING FULL DATA MODE***")
  list_of_seurat_objs <- lapply(names(list_of_matrix_data), Seurat.CreateObj, list_of_matrix_data)
  print(list_of_seurat_objs)
}


#calculate mitochondrial percentage
list_of_seurat_objs <- lapply(list_of_seurat_objs, Seurat.AddMitoRiboPct)

# Filter cells
list_of_seurat_objs <- lapply(list_of_seurat_objs, Seurat.Filter)

# Normalize, scale, find variable genes
list_of_seurat_objs <- lapply(list_of_seurat_objs, Seurat.NormAndScale)
list_of_seurat_objs

#Add name condition in metadata
list_of_seurat_objs[[1]]@meta.data$condition <- condition1
list_of_seurat_objs[[2]]@meta.data$condition <- condition2


# ====== SAVE IMAGE ======
# save.image(file="rsession-1.1_after_data_processing.RData")
# ========================

##########################################################################
############################# CLEAN UP ###################################
##########################################################################
print("Cleaning up variables")

rm(matrix.dat1, matrix.dat2)

##########################################################################
############################# CCA and Merge ##############################
##########################################################################
print("Running CCA")

# Union of highly variable genes
hvg.condition1 <- rownames(x = head(x = list_of_seurat_objs[[1]]@hvg.info, n = 1000))
hvg.condition2 <- rownames(x = head(x = list_of_seurat_objs[[2]]@hvg.info, n = 1000))
hvg.union <- union(x = hvg.condition1, y = hvg.condition2)

#Run Canonical Correlation Analysis
seurat_merge <- RunCCA(object = list_of_seurat_objs[[1]], object2 = list_of_seurat_objs[[2]], genes.use = hvg.union, add.cell.id1=condition1, add.cell.id2=condition2) 
# Without add.cell.id1/2  --> ERROR: Duplicate cell names, please provide 'add.cell.id1' and/or 'add.cell.id2' for unique names
# With add.cell.id: "A bug about the alignment for two-datasets" (Oct 20th): https://github.com/satijalab/seurat/issues/184
# --> need to be on dev branch to fix this error.


##########################################################################
####### Subsetting cells and Seurat objects before alignment  ############
##########################################################################
print("Preparing for alignment")


# Calcuate CalcVarExpRatio
# Before we align the subspaces, we first search for cells whose expression profile cannot be well-explained by low-dimensional CCA, compared to low-dimensional PCA.
seurat_merge <- CalcVarExpRatio(object = seurat_merge, reduction.type = "pca", grouping.var = "condition", dims.use = 1:N_DIMS)

# SPLIT SEURAT OBJECTS 
# seurat_merge.align: Remove Cells which explain too much variance within set (data-set specific cells)
seurat_merge.pre_align <- seurat_merge # save object (copy)
seurat_merge.align <- SubsetData(object = seurat_merge, subset.name = "var.ratio.pca", accept.low = 0.5) # cells for alignment
seurat_merge.discard <- SubsetData(object = seurat_merge, subset.name = "var.ratio.pca", accept.high = 0.5) # discarded cells.


##########################################################################
############################# ALIGNMENT ##################################
##########################################################################
print("Running alignment, tSNE and cluster")

#seurat_merge.align: Align datasets by selected CCA dimensions
seurat_merge.align <- AlignSubspace(object = seurat_merge.align, reduction.type = "cca", grouping.var = "condition", dims.align = 1:N_DIMS)

# Run tSNE and cluster
seurat_merge.align <- RunTSNE(object = seurat_merge.align, reduction.use = "cca.aligned", dims.use = 1:N_DIMS, do.fast = TRUE)
for (tmp.res in res2calculate) {
  seurat_merge.align <- FindClusters(object = seurat_merge.align, reduction.type = "cca.aligned", dims.use = 1:N_DIMS, save.SNN = TRUE, resolution=tmp.res)
}

seurat_merge.align <- SetAllIdent(seurat_merge.align, id=paste0("res.", res.primary)) # *IMPORTANT* set ident to the primary ident. This is needed for correct use of cluster marker gene identification

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", output_prefix))
# ========================

##########################################################################
########################### CLUSTER IDs #################################
##########################################################################
cluster_ids <- sort(unique(as.numeric(seurat_merge.align@meta.data[,paste0("res.",res.primary)]))) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)

##########################################################################
########### Cluster MARKERS: AllMarkers and Conserved Markers ############
##########################################################################
print("Finding markers")

#Identify cluster markers
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.all_markers <- parLapply(cl, cluster_ids, function(x) FindMarkers(object = seurat_merge.align, ident.1 = x, min.pct = 0.25, only.pos = TRUE, thresh.use = 0.25))
##df.cluster_markers <- FindAllMarkers(seurat_merge.align, only.pos=TRUE, min.pct=0.25,  thresh.use=0.25) # TAKES A LONG TIME/PARALLELIZED CODE ABOVE
stopCluster(cl)
list_of_dfs.all_markers <- lapply(list_of_dfs.all_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
names(list_of_dfs.all_markers) <- cluster_ids # set names so bind_rows will get the .id correct
df.cluster_markers <- bind_rows(list_of_dfs.all_markers,.id="cluster") # combine list of dfs into a single data frame

file.out.tmp <- sprintf("%s-all_markers-res_%s.csv", output_prefix, res.primary)
write.csv(df.cluster_markers, file=file.out.tmp, row.names=F, quote=F)


save.image(file=sprintf("%s-rsession_complete.RData", output_prefix))


# Prepare script for parallelization
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()

#Identify Conserved Markers that are present in BOTH conditions (grouping.var="condition")
# we find makers for each cluster: 1 vs rest [ident.2 is NULL (default), then compare to all other cells]
list_of_dfs.conserved_markers <-parLapply(cl, cluster_ids, function(x) FindConservedMarkers(seurat_merge.align, ident.1 = x, grouping.var="condition")) # returns a list of data frames (not named - only 'indexes')
stopCluster(cl)
list_of_dfs.conserved_markers <- lapply(list_of_dfs.conserved_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
names(list_of_dfs.conserved_markers) <- cluster_ids # set names so bind_rows will get the .id correct
df.conserved_markers <- bind_rows(list_of_dfs.conserved_markers,.id="cluster") # combine list of dfs into a single data frame

file.out.tmp <- sprintf("%s-conserved_markers-res_%s.csv", output_prefix, res.primary)
write.csv(df.conserved_markers, file=file.out.tmp ,row.names=F, quote=F)

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", output_prefix))
# ========================

##########################################################################
################## Condition all cells (GLOBAL) DE genes #################
##########################################################################
print("Condition DE genes GLOBAL")

# Identify DE genes across all clusters
df.condition_de_genes_all_cells <- FindMarkers(SetAllIdent(seurat_merge.align, id="condition"), ident.1=condition1, ident.2=condition2, test.use = "wilcox") # Returns data frame


##########################################################################
################## Subset Seurat object into clusters ####################
##########################################################################

# Subset seurat object into clusters
list_of_seurat_objs_per_cluster <- lapply(cluster_ids, function(x) SubsetData(seurat_merge.align, ident.use=x)) # create a list of seurat objects that each contains only cells from a single cluster.
# before setting names() a the list, the names() will be NULL.
names(list_of_seurat_objs_per_cluster) <- cluster_ids # **IMPORTANT** if list is not named, the cluster information COULD BE lost in condition_de_genes() because it may return null.
list_of_seurat_objs_per_cluster <- lapply(list_of_seurat_objs_per_cluster, function(x) SetAllIdent(x, "condition")) # set 'ident' to our condition variable.

##########################################################################
################## Condition DE genes within cluster #####################
##########################################################################
print("Condition DE genes within cluster")

# Identify DE genes within clusters
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
# list_of_dfs.condition_de_genes_per_cluster <- lapply(list_of_seurat_objs_per_cluster, condition_de_genes) # NON PARALLEL VERSION. returns a named list of data frames. Names are the names from the input list.
list_of_dfs.condition_de_genes_per_cluster <-parLapply(cl, list_of_seurat_objs_per_cluster, condition_de_genes) # PARALLEL VERSION. returns a list of data frames (not named - only 'indexes')
stopCluster(cl)
names(list_of_dfs.condition_de_genes_per_cluster) <- cluster_ids # NOT IMPORANT (although lapply will have wrong 'offset' for cluster, i.e. 1, 2, 3 instead of 0, 1, 2, ..): adding names before calling bind_rows. In this way we ensure that the cluster information is kept. this works even if condition_de_genes returns one or more "NULL" values.
df.condition_de_genes_per_cluster <- bind_rows(list_of_dfs.condition_de_genes_per_cluster, .id="cluster") # returns results for all genes (no p_val_adj filter)

file.out.tmp <- sprintf("%s-condition_de_genes_per_cluster-res_%s.csv", output_prefix, res.primary)
write.csv(df.condition_de_genes_per_cluster, file=file.out.tmp ,row.names=F, quote=F)



##########################################################################
########### SAVE FINAL SESSION VARIABLES AND SEURAT OBJECTS  #############
##########################################################################
print("Saving session")

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", output_prefix))
# ========================





















#================================================================================================================================================#
#=========================================================== NOT ALIGNED ANALYSIS  ==============================================================#
#================================================================================================================================================#


print("========== RUNNING Pre-alignment ANALYSIS ==============")

##########################################################################
################### Pre-alignment calculations ###########################
##########################################################################
print("Running pre-alignment calculations: PCA, tSNE, FindClusters")

# seurat_merge.pre_align: Run PCA, tSNE and clustering
# we can use this to illustrate the effect of alignment later
# we run clustering and tsne on PCA space (and not CCA space), because we want to compare the 'standard method' to the 'alignment method'.
seurat_merge.pre_align <- RunPCA(object = seurat_merge.pre_align, pcs.compute = 40, do.print=F) # maxit = 500, weight.by.var = FALSE
seurat_merge.pre_align <- RunTSNE(object = seurat_merge.pre_align, reduction.use = "pca", dims.use = 1:N_DIMS, do.fast = TRUE)
for (tmp.res in res2calculate) {
  seurat_merge.pre_align <- FindClusters(object = seurat_merge.pre_align, reduction.type = "pca", dims.use = 1:N_DIMS, save.SNN = TRUE, resolution=tmp.res)
}
seurat_merge.pre_align <- SetAllIdent(seurat_merge.pre_align, id=paste0("res.", res.primary)) # # set ident to the primary ident

##########################################################################
########################### CLUSTER IDs #################################
##########################################################################
cluster_ids.pre_align <- sort(unique(as.numeric(seurat_merge.pre_align@meta.data[,paste0("res.",res.primary)]))) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)


##########################################################################
################## Condition all cells (GLOBAL) DE genes #################
##########################################################################
print("Condition DE genes GLOBAL")

# Identify DE genes across all clusters
df.condition_de_genes_all_cells.pre_align <- FindMarkers(SetAllIdent(seurat_merge.pre_align, id="condition"), ident.1=condition1, ident.2=condition2, test.use = "wilcox") # Returns data frame

##########################################################################
################## Subset Seurat object into clusters ####################
##########################################################################

# Subset seurat object into clusters
list_of_seurat_objs_per_cluster.pre_align <- lapply(cluster_ids.pre_align, function(x) SubsetData(seurat_merge.pre_align, ident.use=x)) # create a list of seurat objects that each contains only cells from a single cluster.
names(list_of_seurat_objs_per_cluster.pre_align) <- cluster_ids.pre_align # **IMPORTANT** if list is not named, the cluster information COULD BE lost in condition_de_genes() because it may return null.
list_of_seurat_objs_per_cluster.pre_align <- lapply(list_of_seurat_objs_per_cluster.pre_align, function(x) SetAllIdent(x, "condition")) # set 'ident' to our condition variable.


##########################################################################
################## Condition within cluster DE genes  ###################
##########################################################################
print("Condition DE genes within cluster")

# # Identify DE genes within clusters | NON-PARALLEL
# list_of_dfs.condition_de_genes_per_cluster.pre_align <- lapply(list_of_seurat_objs_per_cluster.pre_align, condition_de_genes) # returns a named list of data frames. Names are the names from the input list.
# df.condition_de_genes_per_cluster.pre_align <- bind_rows(list_of_dfs.condition_de_genes_per_cluster.pre_align, .id="cluster")
# df.condition_de_genes_per_cluster.pre_align$cluster <- as.numeric(df.condition_de_genes_per_cluster.pre_align$cluster) - 1 # CORRECTING FOR OFF-SET (Seurat cluster IDs are zero-based. ONLY RUN THIS LINE ONCE

# Identify DE genes within clusters
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.condition_de_genes_per_cluster.pre_align <-parLapply(cl, list_of_seurat_objs_per_cluster.pre_align, condition_de_genes) # PARALLEL VERSION. returns a list of data frames (not named - only 'indexes')
stopCluster(cl)
names(list_of_dfs.condition_de_genes_per_cluster.pre_align) <- cluster_ids.pre_align # NOT IMPORANT (although lapply will have wrong 'offset' for cluster, i.e. 1, 2, 3 instead of 0, 1, 2, ..): adding names before calling bind_rows. In this way we ensure that the cluster information is kept. this works even if condition_de_genes returns one or more "NULL" values.
df.condition_de_genes_per_cluster.pre_align <- bind_rows(list_of_dfs.condition_de_genes_per_cluster.pre_align, .id="cluster") # returns results for all genes (no p_val_adj filter)


file.out.tmp <- sprintf("%s-condition_de_genes_per_cluster_pre_align-res_%s.csv", output_prefix, res.primary)
write.csv(df.condition_de_genes_per_cluster.pre_align, file=file.out.tmp ,row.names=F, quote=F)

##########################################################################
########### SAVE FINAL SESSION VARIABLES AND SEURAT OBJECTS  #############
##########################################################################
print("Saving session - WITH PREALIGNMENT")

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", output_prefix))
# ========================



#========================================================================#
#============================== FINISH ==================================#
#========================================================================#

print("Script DONE!")
