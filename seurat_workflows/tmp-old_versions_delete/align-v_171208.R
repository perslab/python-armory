
# time Rscript alignfinal-v_master.R |& tee alignfinal-v_master.out.txt
  # Bash version > 4

setwd("/home/cbmr/djw472/ygg-projects/timshel/sc-arc_lira/src/align_hfd_ad_lib_lira-vs-hfd_pair_fed-master") # let's not save the 'wd' variable

######################################################################


library(Seurat)
library(dplyr)

library(parallel)


######################################################################
########################### CONSTANTS ################################
######################################################################


# analysis_prefix will be used as prefix for output files.
analysis_prefix <- "lira_vs_pf" 
# analysis_prefix <- "hfd_vs_chow"
# analysis_prefix <- "hfd_vs_pf"

condition1 <-"hfd_ad_lib_lira"
condition2 <-"hfd_pair_fed"

# condition1 <-"hfd_ad_lib"
# condition2 <-"chow"

# condition1 <-"hfd_ad_lib"
# condition2 <-"hfd_pair_fed"



res.primary <- 0.8 # the primary resolution - will be used for cluster marker identification.
res2calculate <- c(0.4, 0.6, 0.8, 1.2, 1.6, 2, 4, 10) # resolutions to calculate for FindClusters
  # OBS: res.primary must be included in res2calculate
  # "name" in metadata is set via paste0("res.", r), so res=2.0 will be "res.2" and NOT "res.2.0"
  # sprintf("res.%s", 2.0) --> "res.2"
  # paste0("res.", 2.0) --> "res.2"

N_DIMS <- 20 # dims.use for CCA (and tSNE)
N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS <- 20 # used for condition_de_genes(seurat_obj)
N_CORES <- 3


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

# Create seurat objects for each experiment
list_of_seurat_objs <- lapply(names(list_of_matrix_data), Seurat.CreateObj, list_of_matrix_data)

# Create seurat objects for each experiment | TESTING
# list_of_seurat_objs <- lapply(names(list_of_matrix_data), Seurat.CreateObj.test_small_dataset, list_of_matrix_data) # TEST on small data set
# list_of_seurat_objs

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
save.image(file=sprintf("%s-rsession_complete.RData", analysis_prefix))
# ========================

##########################################################################
########################### CLUSTER IDs #################################
##########################################################################
cluster_ids <- sort(unique(as.numeric(seurat_merge.align@meta.data[,paste0("res.",res.primary)]))) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)

##########################################################################
########### Cluster MARKERS: AllMarkers and Conserved Markers ############
##########################################################################
print("Finding markers")

#Identify Markers
df.cluster_markers <- FindAllMarkers(seurat_merge.align, only.pos=TRUE, min.pct=0.25,  thresh.use=0.25) # TAKES A LONG TIME
file.out.tmp <- sprintf("%s-all_markers-res_%s.csv", analysis_prefix, res.primary)
write.csv(df.cluster_markers, file=file.out.tmp, row.names=F, quote=F)

save.image(file=sprintf("%s-rsession_complete.RData", analysis_prefix))

# Prepare script for parallelization
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))

#Identify Conserved Markers that are present in BOTH conditions (grouping.var="condition")
# we find makers for each cluster: 1 vs rest [ident.2 is NULL (default), then compare to all other cells]
list_of_dfs.conserved_markers <-parLapply(cl, cluster_ids, function(x) FindConservedMarkers(seurat_merge.align, ident.1 = x, grouping.var="condition")) # returns a list of data frames (not named - only 'indexes')
stopCluster(cl)
list_of_dfs.conserved_markers <- lapply(list_of_dfs.conserved_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
df.conserved_markers <- bind_rows(list_of_dfs.conserved_markers,.id="cluster_no") # combine list of dfs into a single data frame
df.conserved_markers$cluster_no <- as.numeric(df.conserved_markers$cluster_no) - 1 # CORRECTING FOR OFF-SET (Seurat cluster IDs are zero-based. ONLY RUN THIS LINE ONCE

file.out.tmp <- sprintf("%s-conserved_markers-res_%s.csv", analysis_prefix, res.primary)
write.csv(df.conserved_markers, file=file.out.tmp ,row.names=F, quote=F)

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", analysis_prefix))
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
names(list_of_seurat_objs_per_cluster) <- cluster_ids # **IMPORTANT** if list is not named, the cluster information COULD BE lost in condition_de_genes() because it may return null.
list_of_seurat_objs_per_cluster <- lapply(list_of_seurat_objs_per_cluster, function(x) SetAllIdent(x, "condition")) # set 'ident' to our condition variable.


##########################################################################
################## Condition DE genes within cluster #####################
##########################################################################
print("Condition DE genes within cluster")

condition_de_genes <- function(seurat_obj) {
  # Find differentially eseurat_objpressed genes between the two conditions
  # *OBS* note that this function will return NULL if the 'if statement' is not true.
  
  if(sum(seurat_obj@meta.data$condition == condition1) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS 
     & sum(seurat_obj@meta.data$condition == condition2) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS) {
    FindMarkers(seurat_obj, condition1, condition2, test.use = "wilcox") # Returns data frame
      # OBS: SetAllIdent(seurat_obj, "condition") MUST have been called on the Seurat object prior to this, or you will get an error "Identity : hfd_ad_lib_lira not found."
  }
}

# Identify DE genes within clusters
list_of_dfs.condition_de_genes_per_cluster <- lapply(list_of_seurat_objs_per_cluster, condition_de_genes) # returns a named list of data frames. Names are the names from the input list.
list_of_dfs.condition_de_genes_per_cluster <- lapply(list_of_dfs.condition_de_genes_per_cluster, function(x) cbind(x,gene=rownames(x))) # add gene name as column
df.condition_de_genes_per_cluster <- bind_rows(list_of_dfs.condition_de_genes_per_cluster, .id="cluster_no") # returns results for all genes (no p_val_adj filter)
df.condition_de_genes_per_cluster$cluster_no <- as.numeric(df.condition_de_genes_per_cluster$cluster_no) - 1 # CORRECTING FOR OFF-SET (Seurat cluster IDs are zero-based. ONLY RUN THIS LINE ONCE

file.out.tmp <- sprintf("%s-condition_de_genes_per_cluster-res_%s.csv", analysis_prefix, res.primary)
write.csv(df.condition_de_genes_per_cluster %>% filter(p_val_adj < 0.1), file=file.out.tmp ,row.names=F, quote=F)


##########################################################################
########### SAVE FINAL SESSION VARIABLES AND SEURAT OBJECTS  #############
##########################################################################
print("Saving session")

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", analysis_prefix))
# ========================


#========================================================================#
#========================== NOT ALIGNED ANALYSIS  =======================#
#========================================================================#
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

# Identify DE genes within clusters
list_of_dfs.condition_de_genes_per_cluster.pre_align <- lapply(list_of_seurat_objs_per_cluster.pre_align, condition_de_genes) # returns a named list of data frames. Names are the names from the input list.
list_of_dfs.condition_de_genes_per_cluster.pre_align <- lapply(list_of_dfs.condition_de_genes_per_cluster.pre_align, function(x) cbind(x,gene=rownames(x))) # add gene name as column
df.condition_de_genes_per_cluster.pre_align <- bind_rows(list_of_dfs.condition_de_genes_per_cluster.pre_align, .id="cluster_no")
df.condition_de_genes_per_cluster.pre_align$cluster_no <- as.numeric(df.condition_de_genes_per_cluster.pre_align$cluster_no) - 1 # CORRECTING FOR OFF-SET (Seurat cluster IDs are zero-based. ONLY RUN THIS LINE ONCE

##########################################################################
########### SAVE FINAL SESSION VARIABLES AND SEURAT OBJECTS  #############
##########################################################################
print("Saving session - WITH PREALIGNMENT")

# ====== SAVE IMAGE ======
save.image(file=sprintf("%s-rsession_complete.RData", analysis_prefix))
# ========================


#========================================================================#
#============================== FINISH ==================================#
#========================================================================#

print("Script DONE!")
