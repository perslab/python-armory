
library(Matrix)
library(Seurat)
library(dplyr)
library(parallel)

# Rscript <RSCRIPT_FILENAME>.R | tee <RSCRIPT_FILENAME>.out.txt

######################################################################

### DO NOT SET WD FOR THE SCRIPT
# wd <- "/home/cbmr/djw472/ygg-projects/timshel/sc-arc_lira/src"
# setwd(wd)

##set variables
#change names based on files loaded
cond1 <-"hfd_ad_lib_lira"
cond2 <-"hfd_pair_fed"
res <- 0.8
N_DIMS=20 # dims.use for CCA (and tSNE)
N_CORES <- 10

######################################################################


##Useful Functions

#clean up and calculate mito data for each experiment
mitop<-function(x){
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = x@data), ignore.case = TRUE, value = TRUE)
  percent.mito <- Matrix::colSums(x@raw.data[mito.genes, ])/Matrix::colSums(x@raw.data)
  x <- AddMetaData(object = x, metadata = percent.mito, col.name = "percent.mito")
}


#Normalize,scale,find var genes
finalstep<-function(x){
  x<-NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4)
  x<-ScaleData(x, vars.to.regress = c("nUMI","percent.mito"))
  #x<-ScaleData(x, vars.to.regress = c("nUMI","percent.mito"), use.umi=T, model.use="negbinom") # too slow
  x<-FindVariableGenes(x)
}


dat1 <- Read10X(paste0("/raid5/data/sc-10x/data-runs/170612-perslab-arc_lira/",sprintf("agg-arc_lira-%s-no_norm", cond1),"/outs/filtered_gene_bc_matrices_mex/mm10/"))
dat2 <- Read10X(paste0("/raid5/data/sc-10x/data-runs/170612-perslab-arc_lira/",sprintf("agg-arc_lira-%s-no_norm", cond2),"/outs/filtered_gene_bc_matrices_mex/mm10/"))

#data<-list(fgf1=dat1,pf=dat2)
data<-list(cond1=dat1,cond2=dat2)

#create seurat objects for each experiment
seuratobj<-sapply(names(sapply(data, names)), function(x) CreateSeuratObject(raw.data=data[[x]],min.cells=3,min.genes=500, project=x))

#calculate mitochondrial percentage
seuratmito<-lapply(seuratobj, function(x)  mitop(x))

#make plots to determine filtering conditions
sapply(names(sapply(seuratmito, names)), function(x){
  pdf(file=paste0(x,"filter.pdf"))
  par(mfrow=c(1,2))
  GenePlot(object = seuratmito[[x]], gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = seuratmito[[x]], gene1 = "nUMI", gene2 = "nGene")
  dev.off()
})

#filter cells based on plots
seuratfilters <- lapply(seuratmito, function(x) FilterCells(object = x, subset.names = c("nGene", "percent.mito"), 
                                                            low.thresholds = c(200, -Inf), high.thresholds = c(4000,.3)))
# Prepare script for parallelization
#no_cores <- 10
#cl <- makeCluster(no_cores, type = "FORK")
#clusterExport(cl, "finalstep")
#clusterEvalQ(cl, library(Seurat))

#normalize, scale, find variable genes
seurat_finalstep <-lapply(seuratfilters,function(x)finalstep(x))
#stopCluster(cl)

# ====== SAVE IMAGE ======
save.image(file="rsession-1.1_after_data_processing.RData")
# ========================


#Set Variables
dat1<-seurat_finalstep[[1]]
dat2<-seurat_finalstep[[2]]

hvg.cond1 <- rownames(x = head(x = dat1@hvg.info, n = 1000))
hvg.cond2 <- rownames(x = head(x = dat2@hvg.info, n = 1000))
hvg.union <- union(x = hvg.cond1, y = hvg.cond2)

#Rename treatment in metadata
dat1@meta.data$treatment<-cond1
dat2@meta.data$treatment<-cond2


#Run Canonical Correlation Analysis
#seurat_merge <- RunCCA(object = dat1, object2 = dat2, genes.use = hvg.union)
seurat_merge <- RunCCA(object = dat1, object2 = dat2, genes.use = hvg.union, add.cell.id1=cond1, add.cell.id2=cond2) 
# Without add.cell.id1/2  --> ERROR: Duplicate cell names, please provide 'add.cell.id1' and/or 'add.cell.id2' for unique names
# With add.cell.id: "A bug about the alignment for two-datasets" (Oct 20th): https://github.com/satijalab/seurat/issues/184
# --> need to be on dev branch to fix this error.
pdf(paste0(format(Sys.time(),"%Y-%m-%d"), cond1,"_",cond2,"_seurat_merge_PRECCA.pdf"))
p1 <- DimPlot(object = seurat_merge, reduction.use = "cca", group.by = "treatment", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = seurat_merge, features.plot = "CC1", group.by = "treatment", do.return = TRUE)
plot_grid(p1, p2)
dev.off()


#View genes that separate CCA dimensions
PrintDim(object = seurat_merge, reduction.type = "cca", dims.print = 1:2, genes.print = 10)


#Identify how many CCA dimensions to use
pdf(paste0(format(Sys.time(),"%Y-%m-%d"), "_seurat_merge_", cond1,"_",cond2,"_heatmap1.pdf"))
DimHeatmap(object = seurat_merge, reduction.type = "cca", cells.use = 500, dim.use = 1:9,
           do.balanced = TRUE)
dev.off()
pdf(paste0(format(Sys.time(),"%Y-%m-%d"), "_seurat_merge_", cond1,"_",cond2,"_heatmap2.pdf"))
DimHeatmap(object = seurat_merge, reduction.type = "cca", cells.use = 500, dim.use = 10:20,
           do.balanced = TRUE)
dev.off()

# ====== SAVE IMAGE ======
save.image(file="rsession-1.2_after_CCA.RData")
# ========================


#Remove Cells which explain too much variance within set (data-set specific cells)
seurat_merge <- CalcVarExpRatio(object = seurat_merge, reduction.type = "pca", grouping.var = "treatment",
                                dims.use = 1:N_DIMS)
seurat_merge.save <- seurat_merge
seurat_merge <- SubsetData(object = seurat_merge, subset.name = "var.ratio.pca", accept.low = 0.5)
seurat_merge.discard <- SubsetData(object = seurat_merge, subset.name = "var.ratio.pca", accept.high = 0.5)
save(seurat_merge.discard, file = "dataset_specific_cells")

#Align datasets by selected CCA dimensions
seurat_merge <- AlignSubspace(object = seurat_merge, reduction.type = "cca", grouping.var = "treatment",
                              dims.align = 1:N_DIMS)
p1 <- VlnPlot(object = seurat_merge, features.plot = "ACC1", group.by = "treatment",
              do.return = TRUE)
p2 <- VlnPlot(object = seurat_merge, features.plot = "ACC2", group.by = "treatment",
              do.return = TRUE)
pdf(paste0(format(Sys.time(),"%Y-%m-%d"), "_seurat_merge_", cond1,"_",cond2,"_CCA_alignment.pdf"))
plot_grid(p1, p2)
dev.off()


#Cluster
seurat_merge <- RunTSNE(object = seurat_merge, reduction.use = "cca.aligned", dims.use = 1:N_DIMS,
                        do.fast = TRUE)
seurat_merge <- FindClusters(object = seurat_merge, reduction.type = "cca.aligned", dims.use = 1:N_DIMS,
                             save.SNN = TRUE, resolution=res)
p1 <- TSNEPlot(object = seurat_merge, group.by = "treatment", do.return = TRUE, pt.size = 1)
p2 <- TSNEPlot(object = seurat_merge, do.return = TRUE, pt.size = 1, do.label=TRUE, label.size = 5)

pdf(paste0(format(Sys.time(),"%Y-%m-%d"), "_seurat_", cond1,"_",cond2,"_merge_tsne_",res,".pdf"), width=10,height=7.5)
plot_grid(p1, p2)
dev.off()
save(seurat_merge, file=paste0("seurat_",cond1,"_",cond2,"_merged_res",res))

# ====== SAVE IMAGE ======
save.image(file="rsession-3.1_after_findclusters.RData")
# ========================


# Subset seurat object into clusters
col<-paste0("res.",res)
x<-seq_len(max(as.numeric(seurat_merge@meta.data[,col]))+1)-1
clusters<-sapply(x, function(x) SubsetData(seurat_merge, ident.use=x))
names(clusters)<-paste0("cluster_",x)
clusters<-lapply(clusters, function(x) SetAllIdent(x, "treatment")) # set 'ident' to our treatment variable.
#clusters<-lapply(clusters, function(x) SetAllIdent(x, "orig.ident"))


# Identify DE genes within clusters
N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS <- 20
cluster.markers<-lapply(clusters, function(x) if(sum(x@meta.data$treatment == cond1) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS 
                                                 & sum(x@meta.data$treatment == cond2) > N_MIN_CELLS_PER_TREATMENT_FOR_DE_ANALYSIS) {
  FindMarkers(x, cond1, cond2)})
cluster.markers.de<-lapply(cluster.markers, function(x) cbind(x,gene=rownames(x)))
cluster.markers.de<-bind_rows(cluster.markers.de,.id="id")
decluster<-cluster.markers.de[cluster.markers.de$p_val_adj<.05,]
write.csv(decluster,paste0(cond1,"_",cond2,"_declusters_res",res,".csv"), row.names=F,quote=F)

# ====== SAVE IMAGE ======
save.image(file="rsession-3.2_after_find_DE_within_clusters.RData")
# ========================

# Prepare script for parallelization
cl <- makeCluster(N_CORES, type = "FORK")
clusterEvalQ(cl, library(Seurat))

#Identify Conserved Markers
conservedmarks<-parLapply(cl, x, function(x) FindConservedMarkers(seurat_merge, ident.1 = x, grouping.var="treatment"))
stopCluster(cl)
consmark.top<-lapply(conservedmarks, function(x) cbind(x,gene=rownames(x)))
consmark.top<-bind_rows(consmark.top,.id="id")
write.csv(consmark.top, paste0(cond1,"_",cond2,"_consmark_res",res,".csv"),row.names=F, quote = F)

# ====== SAVE IMAGE ======
save.image(file="rsession-4.1_after_conservedmarkers.RData")
# ========================

#Identify Markers
all_markers <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
write.csv(all_markers,paste0(cond1,"_",cond2,"_allmarkers_res",res,".csv"),row.names=F, quote = F)

# ====== SAVE IMAGE ======
save.image(file="rsession-4.2_after_findallmarkers.RData")
# ========================



