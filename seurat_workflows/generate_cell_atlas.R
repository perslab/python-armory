# Title: Script generate a cell atlas. Will do all the heavy computations that can afterwards be analyzed/plotted in R Notebooks.
# Author: Pascal N. Timshel, Pers Lab
# Date: December 2017

######################################################################
############################## USAGE #################################
######################################################################

# time Rscript generate_cell_atlas.R |& tee generate_cell_atlas.out.txt

######################################################################
############################## SETUP #################################
######################################################################

library(Seurat)
library(tidyverse)
library(Matrix)

library(parallel)


dir.project_src <- "/projects/timshel/sc-arc_lira/src"
source(sprintf("%s/constants-sample_meta_data.R", dir.project_src)) # loads sample meta-data


######################################################################
########################### CONSTANTS ################################
######################################################################

output_prefix <- "arc_lira_cell_atlas"

res.primary <- 0.8 # the primary resolution - will be used for cluster marker identification.
res2calculate <- c(0.4, 0.6, 0.8, 1.2, 1.6, 2, 4, 10) # resolutions to calculate for FindClusters. OBS: res.primary must be included in res2calculate

N_DIMS <- 10 # dims.use for CLUSTERING and tSNE
N_CORES <- 10

dir.10x_data <- "/data/sc-10x/data-runs/170612-perslab-arc_lira/agg-arc_lira-s1_12-no_norm/outs/filtered_gene_bc_matrices_mex/mm10" # path to data

######################################################################
########################### LOAD DATA ################################
######################################################################
print("Loading data...")

data.matrix <- Read10X(data.dir = dir.10x_data)
dim(data.matrix) # genes are rows; cells are columns

seurat_obj <- CreateSeuratObject(data.matrix, min.cells = 3, min.genes = 200, is.expr = 0)


######################################################################
########################### META DATA ################################
######################################################################

### Add metadata
# We use the data in constants-sample_meta_data.R to set the correct annotations.
sample_agg_idx <- as.numeric(sapply(strsplit(seurat_obj@cell.names, split = "-"), '[[', 2)) 
# unique(sample_agg_idx) # --> 1,2,3,...,12
df.metadata <- data.frame(sample_name=META.sample_name[sample_agg_idx],
                          group_name=META.group_name[sample_agg_idx],
                          batch=META.batch[sample_agg_idx], 
                          row.names=seurat_obj@cell.names)

seurat_obj <- AddMetaData(seurat_obj, df.metadata)

# Create metadata columns for annotations and subannotations
seurat_obj@meta.data[,'annotation'] <- NA
seurat_obj@meta.data[,'subannotation'] <- NA




######################################################################
######################## DATA PROCESSING #############################
######################################################################
print("Processing data...")

# Calculate percent ribosomal and mitochondrial genes.
ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = seurat_obj@data), value = TRUE)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_obj@data), value = TRUE, ignore.case=TRUE)
percent.ribo <- Matrix::colSums(seurat_obj@raw.data[ribo.genes, ])/Matrix::colSums(seurat_obj@raw.data)
percent.mito <- Matrix::colSums(seurat_obj@raw.data[mito.genes, ])/Matrix::colSums(seurat_obj@raw.data)
seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.ribo, col.name = "percent.ribo")
seurat_obj <- AddMetaData(object = seurat_obj, metadata = percent.mito, col.name = "percent.mito")


### Filter out cells with few reads and few genes.
# unique gene counts: <200 or >5000
# percent mitochrondria: >25%
seurat_obj <- FilterCells(object=seurat_obj, subset.names=c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.25))

### Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

### Scale data and regress out `nUMI` and `percent.mito`
#seurat_obj <- ScaleData(seurat_obj, model.use = "negbinom", use.umi = TRUE, vars.to.regress = c("nUMI", "percent.mito")) # UMI model is slow...
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nUMI", "percent.mito"))
# TODO: consider regressing out "percent.ribo" and "Rn45s"

######################################################################
############################### PCA ##################################
######################################################################
print("Running PCA...")

### # Principal component analysis
# Find variable genes
seurat_obj <- FindVariableGenes(seurat_obj, do.plot=F)
# Run Principal Component Analysis.
seurat_obj <- RunPCA(object = seurat_obj, pcs.compute = 40, do.print=F) # maxit = 500 , weight.by.var = FALSE,


######################################################################
########################## CLUSTERING ################################
######################################################################
print("Running Clustering...")

## # Cluster cells
for (tmp.res in res2calculate) {
  seurat_obj <- FindClusters(object = seurat_obj, reduction.type = "pca", dims.use = 1:N_DIMS, save.SNN = TRUE, print.output = 0, resolution=tmp.res)
}

seurat_obj <- SetAllIdent(seurat_obj, id=paste0("res.", res.primary)) # *IMPORTANT* set ident to the primary ident. This is needed for correct use of cluster marker gene identification


######################################################################
############################## TSNE ##################################
######################################################################
print("Running tSNE...")

### Run tSNE algorithm to get a cell-type atlas
seurat_obj <- RunTSNE(object = seurat_obj, dims.use = 1:N_DIMS) # perplexity=30, dim.embed = 2


##########################################################################
########################### CLUSTER IDs #################################
##########################################################################
cluster_ids <- sort(unique(as.numeric(seurat_obj@meta.data[,paste0("res.",res.primary)]))) # numeric vector, 0, 1, 2, 3, ... (N_CLUSTERS-1)


######################################################################
######################### MARKER GENES ###############################
######################################################################
print("Finding markers...")

#Identify cluster markers
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(rlang)) # needed for using UQ()
list_of_dfs.all_markers <- parLapply(cl, cluster_ids, function(x) FindMarkers(object = seurat_obj, ident.1 = x, min.pct = 0.25, only.pos = TRUE, thresh.use = 0.25))
##df.cluster_markers <- FindAllMarkers(seurat_merge.align, only.pos=TRUE, min.pct=0.25,  thresh.use=0.25) # TAKES A LONG TIME/PARALLELIZED CODE ABOVE
stopCluster(cl)
list_of_dfs.all_markers <- lapply(list_of_dfs.all_markers, function(x) cbind(x,gene=rownames(x))) # add gene name as column
names(list_of_dfs.all_markers) <- cluster_ids # set names so bind_rows will get the .id correct
df.cluster_markers <- bind_rows(list_of_dfs.all_markers,.id="cluster") # combine list of dfs into a single data frame



##########################################################################
########### SAVE FINAL SESSION VARIABLES AND SEURAT OBJECTS  #############
##########################################################################

print("Saving session...")
save.image(file=sprintf("%s-rsession_complete.RData", output_prefix)) # save image

print("Saving Seurat object...")
save(seurat_obj, df.cluster_markers, file=sprintf("%s-core_objs.RData", output_prefix))

print("DONE with script!")


