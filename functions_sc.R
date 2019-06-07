#' @title sc-analysis functions
#' @author Jonatan Thompson, Pers lab, rkm916 at ku dot dk
#' 
######################################################################
################### FEATUREPLOT WITH 3+ FEATURES #####################
######################################################################

FeaturePlotPlus <- function(seurat_obj, 
                            genes=NULL,
                            metadata.feats=NULL,
                            detection.threshold=1,
                            colors.indiv=NULL,
                            colors.comb=NULL,
                            combinations.exclude=NULL,
                            plot.title,
                            ...) {
  
  #' @usage: wrapper for Seurat::DimPlot that works like FeaturePlot but for several genes
  #' @param seurat_obj: a seurat object with raw.data slot
  #' @param genes: character vector of gene names
  #' @param metadata.feats: list of vectors of metadata levels named by the feature
  #' @param detection.threshold: minimum number of (raw) transcripts for coloring a cell, defaults to 1
  #' @param colors.indiv: List of vectors named 'genes' and 'metadata.feats', giving colors for cells expressing only that feature. If only plotting genes or metadata,
  #' a vector is sufficient. Cells expression a combination of features will show a gradient, if NULL uses default colors.
  #' @param colors.comb: Named list with a color for each combination, e.g. c(Th="blue", Lepr="yellow", Th_Lepr="green")
  #' @param combinations exclude: a list of character vectors of combinations of features to exclude
  #' @param plot.title: passed via TSNEPlot on to ggplot
  #' @param ... : additional named arguments to pass to TSNEPlot onto DimPlot onto ggplot
  #' @value: a ggplot object
  #' @packages: Seurat (including dependencies Matrix and ggplot2), magrittr/dplyr, colorspace
  
  # TODO: Add support for metadata features. Binary over / under threshold?
  # TODO: Add a 'graded' argument, expression valued in 5 bins (and change color scheme accordingly)
  # TODO: Keep track of (normalized) counts and set alpha accordingly when merging (make it an option, though)
  # TODO: Allow the user to pass arguments to DimPlot (see how other wrapper functions do it)
  
  require(Seurat)
  require(magrittr)
  require(colorspace)
  
  if (length(genes) + length(unlist(metadata.feats)) >5) warning("Plotting more than five features is likely to be illegible")
  if (length(genes) + length(unlist(metadata.feats)) >17) stop("Too many features (max 17)")
  if (!is.null(genes)) if (any(!genes %in% rownames(seurat_obj@raw.data))) {
    stop(paste0(paste0(genes[!genes %in% rownames(seurat_obj@raw.data)], collapse= " "), " not found in seurat object raw data"))
  }
  
  #stopifnot(xor(is.null(colors.indiv), is.null(colors.comb)))
  # Get the idx of cells expressing the genes
  #list_idx_gene <- sapply(genes, function(gene) seurat_obj@raw.data[gene,]>= detection.threshold)
  
  # Get all combinations of the genes
  genes_tmp <- c(genes, character(length(genes)-1))
  list_comb <- combn(simplify = F, x = genes_tmp, m=length(genes)) %>% unique
  list_idx_notempty <- lapply(list_comb, function(comb) nchar(comb)>0)
  list_comb <- mapply(function(idx_notempty, comb) comb[idx_notempty], comb=list_comb, idx=list_idx_notempty, SIMPLIFY=F)
  
  # Filter out unwanted combinations
  if (!is.null(combinations.exclude)) list_comb <- Filter(f=function(comb) {
    all(sapply(combinations.exclude, function(comb_exclude) {!all(comb %in% comb_exclude) | !all(comb_exclude %in% comb)}))
  }, list_comb)
  
  names_comb <- sapply(X=list_comb, FUN=function(comb)paste0(comb,collapse="_"))
  
  # test argument
  if (!is.null(colors.comb)) stopifnot(length(names_comb)==length(colors.comb))
  
  names(list_comb) <- names_comb
  names_comb <- sort(names_comb)
  list_comb <- list_comb[order(names(list_comb))]
  
  list_idx_comb <- mapply(function(comb,name) {
    comb_complement <- setdiff(genes, comb)
    
    # Get the intersect of the idx of cells expressing the genes in the combination
    list_idx_gene <- lapply(X=comb, FUN=function(gene) seurat_obj@data[gene, ]>= detection.threshold)
    idx_comb_intersect <- Reduce(x=list_idx_gene, f = '&')
    
    # Get the **union** of the idx of the cells expressing the genes which are **not** in the combination (the complementary set of genes)
    if (length(comb_complement)>= detection.threshold) {
      list_idx_gene_complement <- lapply(X=comb_complement, FUN=function(gene) seurat_obj@data[gene, ]>= detection.threshold)
      idx_comb_complement_union <- Reduce(x=list_idx_gene_complement, f = '|')
    } else {
      idx_comb_complement_union <- logical(length(idx_comb_intersect))
    }
    # Get the idx of cells expressing the combination intersect but not the complementary combination union
    out <- idx_comb_intersect & !idx_comb_complement_union
    names(out) <- name
    out
  }, comb=list_comb, name = names_comb, SIMPLIFY=F)
  
  # Make a single vector of labels 
  vec_labels <- character(length(seurat_obj@ident))
  
  for (i in 1:length(list_idx_comb)) {
    vec_labels[list_idx_comb[[i]]] <- names_comb[i]
  }
  
  vec_labels[nchar(vec_labels)==0] <- "None"#NA_character_
  
  # Add vector of labels to Seurat object as metadata
  df_labels <- data.frame(feats_to_plot = vec_labels, row.names = colnames(seurat_obj@data))
  
  seurat_obj <- AddMetaData(seurat_obj, df_labels)
  
  # Make colors for plotting
  # make a 'primary' palette, supplementing with additional colors if needed
  if (is.null(colors.comb)) {
    
    # Find colors for each combination
    colors_1 <- if (is.null(colors.indiv)) c("tomato", "yellow2", "deepskyblue") else colors.indiv
    
    #if (length(colors_1) != length(genes)) stop("colors.use must have the same length as genes")
    
    if (length(genes)>length(colors_1)) {
      
      colors_avail <- c("antiquewhite3", 
                        "aquamarine2", 
                        "bisque1", 
                        "blueviolet", 
                        "brown3", 
                        "burlywood3", 
                        "coral", 
                        "darkblue", 
                        "darkgoldenrod1", 
                        "darkgreen", 
                        "darkkhaki", 
                        "darkolivegreen4", 
                        "darkorchid2", 
                        "deeppink1") 
      
      colors_constituent <- colors_1
      
      while(length(colors_constituent)<length(genes)) {
        size <- length(genes)-length(colors_constituent)
        colors_extra <- sample(x=colors_avail, size = size,replace = F)
        colors_constituent <- unique(c(colors_1,colors_extra))
      }
      
    } else {
      colors_constituent <- colors_1[1:length(genes)] # start with primary colors, only take rando
    }
    
    names(colors_constituent) <- genes
    
    colors_plot <- vector(mode="list", length=length(list_comb))
    names(colors_plot) <- names_comb
    
    for (i in 1:length(list_comb)) {
      if (length(list_comb[[i]])==1) {
        colors_plot[[i]] <- colors_constituent[match(list_comb[[i]], genes)] %>% col2rgb %>% t %>% RGB
      } else {
        colors_to_combine <- colors_constituent[match(list_comb[[i]], names(colors_constituent))] 
        color_combined <- colors_to_combine[1] %>% col2rgb %>% t %>% RGB
        alpha = 1/(1:length(list_comb[[i]]))
        for (j in 2:length(colors_to_combine)) {
          color2 <- colors_to_combine[j] %>% col2rgb %>% t %>% RGB
          color_combined <- mixcolor(alpha=alpha[j], color1 = color_combined, color2=color2)
        }
        for (k in 1:3) { color_combined@coords[[k]] <- color_combined@coords[[k]]*0.92^length(colors_to_combine)^2}
        colors_plot[[i]] <- color_combined 
      }
      
    }
    
  } else { # if the user has provided a list of colors for each combination
    colors.comb <- colors.comb[match(names_comb,names(colors.comb))]
    colors_plot <- lapply(colors.comb, function(col) {
      col %>% col2rgb %>% t %>% RGB
    })
  }
  # Convert from RGB objects to hex (character)
  colors_plot_hex <- sapply(colors_plot, function(color){
    color_hex <- rgb(red = color@coords[1], green = color@coords[2], blue = color@coords[3], maxColorValue = 255)
    color_hex
  })
  
  # add grey for None
  colors_plot_hex = c(colors_plot_hex, "None"="#e2e9ee")
  
  list_args <- list(...) # take ... arguments
  list_args[["object"]]=seurat_obj
  list_args[["do.label"]] = F
  list_args[["colors.use"]] = colors_plot_hex
  list_args[["do.return"]] = T
  list_args[["group.by"]] = "feats_to_plot"
  list_args[["no.legend"]] = F
  list_args[["plot.title"]] = plot.title
  list_args[["na.value"]] = "grey90"
  
  # Make the plot
  p <- do.call(what = TSNEPlot, 
               args= list_args)
  
  return(p)
}

######################################################################
########################### GENE ORTHOLOGY REMAP #####################
######################################################################

#          from_organism = c("mmusculus", "hsapiens"),
#          to_organism = c("mmusculus", "hsapiens"))
#          orthology_mapping_path="/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz",
#          synonym_mapping_path="/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz",

######################################################################
############################### DOUBLETFINDER ########################
######################################################################

# 20190507 DEPRECATED 
# changes to original:
# * instead of using expected.doublets as cut-off, uses outlier criterion pANN > Q3+1.5*(Q3-Q1)
# * rename artificial doublets "thisisfake" rather than "X" to avoid issues e.g. with 10X cells!
# * fix a bug with PCA

# doubletFinder = function(seu, 
#                          proportion.artificial = 0.25, 
#                          proportion.NN = 0.02) 
# {
#   # if (expected.doublets == 0) {
#   #   stop("Need to set number of expected doublets...")
#   # }
#   print("Creating artificial doublets...")
#   data <- GetAssayData(seu, slot="counts")[, colnames(seu)]
#   real.cells <- colnames(seu)
#   n_real.cells <- length(real.cells)
#   n_doublets <- round(n_real.cells/(1 - proportion.artificial) - 
#                         n_real.cells) # these will be artificial doublets
#   real.cells1 <- sample(real.cells, n_doublets, replace = TRUE) # draw samples of real cells of size of n artificial doublets to make
#   real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
#   doublets <- (data[, real.cells1] + data[, real.cells2])/2 # make artificial doublets
#   colnames(doublets) <- paste("thisisfake", 1:n_doublets, sep = "")
#   data_wdoublets <- cbind(data, doublets)
#   print("Creating Seurat object...")
#   seu_wdoublets <- Seurat::CreateSeuratObject(counts = data_wdoublets)
#   print("Normalizing Seurat object...")
#   seu_wdoublets <- Seurat::NormalizeData(seu_wdoublets)
#                                          #normalization.method = seu@calc.params$NormalizeData$normalization.method, 
#                                          #scale.factor = seu@calc.params$NormalizeData$scale.factor)
#   print("Finding variable genes...")
#   seu_wdoublets <- Seurat::FindVariableFeatures(seu_wdoublets)#, 
#                                              #do.plot = FALSE, 
#                                              #x.low.cutoff = seu@calc.params$FindVariableGenes$x.low.cutoff, 
#                                              #x.high.cutoff = seu@calc.params$FindVariableGenes$x.high.cutoff, 
#                                              #y.high.cutoff = seu@calc.params$FindVariableGenes$y.high.cutoff, 
#                                              #y.cutoff = seu@calc.params$FindVariableGenes$y.cutoff)
#   print("Scaling data...")
#   seu_wdoublets <- Seurat::ScaleData(seu_wdoublets, block.size=15000, min.cells.to.block = 5000)
#   
#   print("Running PCA...")
#   seu_wdoublets <- tryCatch({
#     pcs.compute1 <- min(40, min(ncol(GetAssayData(seu_wdoublets, slot="scale.data")), length(VariableFeatures(seu_wdoublets)))%/%2)
#     Seurat::RunPCA(seu_wdoublets, 
#                    #pc.genes = seu_wdoublets@var.genes, 
#                    do.print = F,
#                    pcs.compute = pcs.compute1,
#                    verbose=T)
#   }, error= function(err) {
#     pcs.compute2 <- min(20, min(ncol(seu_wdoublets@scale.data), length(seu_wdoublets@var.genes))%/%3)
#     warning(paste0("RunPCA failed with ", pcs.compute1," components, trying again with ", pcs.compute2, " components"))
#     Seurat::RunPCA(seu_wdoublets, 
#                    #pc.genes = seu_wdoublets@var.genes, 
#                    do.print = F,
#                    pcs.compute = pcs.compute2,
#                    verbose=T)}) 
#   
#   cell.names <- colnames(seu_wdoublets)
#   nCells <- length(cell.names)
#   print("Calculating PC distance matrix...")
#   #PCs <- 1:min(max(seu@calc.params$RunTSNE$dims.use), ncol(seu_wdoublets@dr$pca@cell.embeddings))
#   PCs <- 1:ncol(Embeddings(object=seu_wdoublets, reduction="pca"))
#   if (length(PCs) == 0) {
#    stop("Need to run tSNE on original Seurat object...")
#   }
#   pca.coord <- Embeddings(object=seu_wdoublets, reduction="pca")#seu_wdoublets@dr$pca@cell.embeddings[, PCs]
#   rm(seu_wdoublets)
#   gc()
#   dist.mat <- as.matrix(dist(pca.coord)) # distances in PCA space
#   dist.mat <- dist.mat[, -grep("thisisfake", colnames(dist.mat))] # so keep fake cell rows but not columns
#   pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
#   rownames(pANN) <- real.cells
#   colnames(pANN) <- "pANN"
#   k <- round(nCells * proportion.NN) # k is how many neighbours to check
#   for (i in 1:n_real.cells) { # i goes up to the last real cell in the distance matrix
#     neighbors <- order(dist.mat[, i]) # all cell numbers ranked by distance (i.e. closest first)
#     neighbors <- neighbors[2:(k + 1)] # get k nearest neighbours
#     neighbor.names <- rownames(dist.mat)[neighbors] # these will be cell names. Artificial cells have "X" in them
#     pANN[i, 1] <- length(grep("thisisfake", neighbor.names))/k # for real cell i, how big a proportion of k nearest neighbours are artificial?
#   }
#   seu[["pANN"]] <- pANN # add props fake neighbours to original (real) seurat obj
#   predictions <- as.data.frame(rep("Singlet", n_real.cells), 
#                                ncol = 1, row.names = real.cells , stringsAsFactors = FALSE) # initialise prediction metadata column
#   ### modified @author Jonatan Thompson, jjt3f2188@gmail.com  @date 181108 ###
#   Q1 <- quantile(x = pANN[,1], probs=0.25)
#   Q3 <- quantile(x = pANN[,1], probs=0.75)
#   pANN_inter_Q_range <- as.numeric(Q3-Q1)
#   doublet.predictions <- colnames(seu)[pANN>Q3+1.5*(Q3-Q1)]
#   # doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pANN, 
#   #                                                      decreasing = TRUE)] # order real cells by prop fake neighbours 
#   # doublet.predictions <- doublet.predictions[1:expected.doublets] # take the top of the list, based on absolute number
#   # 
#   predictions[doublet.predictions, ] <- "Doublet"
#   colnames(predictions) <- "pANNPredictions"
#   seu[["pANNPredictions"]] <-  predictions
#   return(seu)
# }


######################################################################
#################### BACKGROUND FINDER ###############################
######################################################################

# backgroundFinder = function(seu, 
#                          background, 
#                          proportion.NN = 0.02) 
# {
#   # if (expected.doublets == 0) {
#   #   stop("Need to set number of expected doublets...")
#   # }
#   #print("Creating artificial doublets...")
#   data <- seu@raw.data[, seu@cell.names]
#   real.cells <- seu@cell.names
#   n_real.cells <- length(real.cells)
#   #n_doublets <- round(n_real.cells/(1 - proportion.artificial) - 
#   #                      n_real.cells) # these will be artificial doublets
#   #real.cells1 <- sample(real.cells, n_doublets, replace = TRUE) # draw samples of real cells of size of n artificial doublets to make
#   #real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
#   #doublets <- (data[, real.cells1] + data[, real.cells2])/2 # make artificial doublets
#   colnames(background) <- paste("background_",colnames(background), sep = "")
#   data_wbackground <- cbind(data, background)
#   print("Creating Seurat object...")
#   seu_wbackground <- Seurat::CreateSeuratObject(raw.data = data_wbackground)
#   print("Normalizing Seurat object...")
#   seu_wbackground <- Seurat::NormalizeData(seu_wbackground)
#   print("Finding variable genes...")
#   seu_wbackground <- Seurat::FindVariableGenes(seu_wbackground, 
#                                              do.plot = FALSE)
#   print("Scaling data...")
#   seu_wbackground <- Seurat::ScaleData(seu_wbackground, display.progress = TRUE)
#   print("Running PCA...")
#   seu_wbackground <- tryCatch({
#     pcs.compute1 <- min(40, min(ncol(seu_wbackground@scale.data), length(seu_wbackground@var.genes))%/%2)
#     Seurat::RunPCA(seu_wbackground, 
#                    pc.genes = seu_wbackground@var.genes, 
#                    pcs.print = 0,
#                    pcs.compute = pcs.compute1,
#                    verbose=T)
#   }, error= function(err) {
#     pcs.compute2 <- min(20, min(ncol(seu_wbackground@scale.data), length(seu_wbackground@var.genes))%/%3)
#     warning(paste0("RunPCA failed with ", pcs.compute1," components, trying again with ", pcs.compute2, " components"))
#     Seurat::RunPCA(seu_wbackground, 
#                    pc.genes = seu_wbackground@var.genes, 
#                    pcs.print = 0,
#                    pcs.compute = pcs.compute2,
#                    verbose=T)}) 
#   
#   cell.names <- seu_wbackground@cell.names
#   nCells <- length(cell.names)
#   print("Calculating PC distance matrix...")
#   PCs <- 1:ncol(seu_wbackground@dr$pca@cell.embeddings)
#   if (length(PCs) == 0) {
#     stop("Need to run tSNE on original Seurat object...")
#   }
#   pca.coord <- seu_wbackground@dr$pca@cell.embeddings[, PCs]
#   rm(seu_wbackground)
#   gc()
#   dist.mat <- as.matrix(dist(pca.coord)) # distances in PCA space
#   dist.mat <- dist.mat[, -grep("_background", colnames(dist.mat))] # so keep fake cell rows but not columns
#   pbackground <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
#   rownames(pbackground) <- real.cells
#   colnames(pbackground) <- "pbackground"
#   k <- round(nCells * proportion.NN) # k is how many neighbours to check
#   for (i in 1:n_real.cells) { # i goes up to the last real cell in the distance matrix
#     neighbors <- order(dist.mat[, i]) # all cell numbers ranked by distance (i.e. closest first)
#     neighbors <- neighbors[2:(k + 1)] # get k nearest neighbours
#     neighbor.names <- rownames(dist.mat)[neighbors] # these will be cell names. Artificial cells have "X" in them
#     pbackground[i, 1] <- length(grep("_background", neighbor.names))/k # for real cell i, how big a proportion of k nearest neighbours are artificial?
#   }
#   seu <- Seurat::AddMetaData(seu, metadata = pbackground, col.name = "pbackground") # add props fake neighbours to original (real) seurat obj
#   predictions <- as.data.frame(rep("real_cell", n_real.cells), 
#                                ncol = 1, stringsAsFactors = FALSE) # initialise prediction metadata column
#   rownames(predictions) <- real.cells 
#   ### modified @author Jonatan Thompson, jjt3f2188@gmail.com  @date 181108 ###
#   pbackgroundsummary <- summary(seu@meta.data$pbackground) 
#   Q1 <- as.numeric(pbackgroundsummary[2])
#   Q3 <- as.numeric(pbackgroundsummary[5])
#   pbackground_inter_Q_range <- as.numeric(Q3-Q1)
#   background.predictions <- rownames(seu@meta.data)[pbackground>Q3+1.5*(Q3-Q1)]
#   # doublet.predictions <- rownames(seu@meta.data)[order(seu@meta.data$pbackground, 
#   #                                                      decreasing = TRUE)] # order real cells by prop fake neighbours 
#   # doublet.predictions <- doublet.predictions[1:expected.background] # take the top of the list, based on absolute number
#   # 
#   predictions[background.predictions, ] <- "background"
#   colnames(predictions) <- "background_predictions"
#   seu <- Seurat::AddMetaData(seu, metadata = predictions, 
#                              col.name = "background_predictions")
#   return(seu)
# }

# # test backgroundfinder
# if (FALSE) {
# 
#   data_tmp <- Seurat::Read10X(data.dir = "/nfsdata/data/sc-10x/data-runs/180511-perslab-immunometab/242L-5000_cells/outs/raw_gene_bc_matrices/hg19/")
#   n_cells <- 1749
#   background_start <- 6000
#   n_background <- n_cells%/%3
#   data_tmp %>% colSums -> nUMI_sums
#   nUMI_sums %>% rank -> rank_tmp
#   idx_background <- which(rank_tmp <= (length(rank_tmp)-background_start) & rank_tmp > (length(rank_tmp)-(background_start+n_background))) 
#   background= data_tmp[, idx_background]
#   
#   # TODO: need to process bacground befpre
#   
# }

######################################################################
############################ GENE MAP ################################
######################################################################

gene_map <- function(dataIn,
                     colGene = NULL,
                     df_mapping,
                     from="hgnc", 
                     to="ensembl",
                     replace = F,
                     na.rm = T) {
  #' @usage map genes in a dataframe column, vector, or list of vectors between naming schemes in df_mapping
  #' @param dataIn data.frame with a column or rownames containing genes to remap,
  #' matrix with gene rownames, a list of vectors, or a vector, either named numeric or character 
  #' If a list, if the vectors are numeric, the vector names are assumed to be genes, if the vectors are character
  #' the vector values are assumed to be genes
  #' @param colGene if dataIn is a dataframe, the name of the gene column; in this case NULL implies rownames 
  #' @param df_mapping data.frame or matrix with columns corresponding to 'from' and 'to' arguments (using grep partial matching)
  #' @param from df_mapping colnames
  #' @param to df_mapping colnames
  #' @param replace boolean; if dataIn is a data.frame, TRUE replaces original gene names, FALSE adds a new column to the data.frame
  #' @param na.rm boolean; remove genes that fail to map or leave them as NAs?; defaults to TRUE
  #' @value an object of the same format as dataIn with new gene names
  
  #' @example : 
  # df_test <- gene_map(dataIn=load_obj("/projects/jonatan/tmp-mousebrain/tables/mousebrain_Vascular_ClusterName_1_PER3_kIM.csv"),
  #                     colGene ="genes",
  #                     df_mapping=load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"),
  #                     from="ensembl", 
  #                     to="gene_name_optimal",
  #                     replace=T,
  #                     na.rm=T)
  # orthologue mapping:
  # /projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz
  
  stopifnot(class(df_mapping)=="data.frame")
  stopifnot(length(from)>0 & length(to)>0)
  
  fromMapCol <- if (from %in% colnames(df_mapping)) from else grep(pattern=from, x = colnames(df_mapping), ignore.case=T, value = T)
  toMapCol <- if (to %in% colnames(df_mapping)) to else grep(pattern=to, x = colnames(df_mapping), ignore.case=T, value = T)
  
  if (length(fromMapCol)==0) stop(paste0(from, " not found in df_mapping column names"))
  if (length(toMapCol)==0) stop(paste0(to, " not found in df_mapping column names"))

  if (!is.null(dim(dataIn))) {
    
    if(is.null(colGene) & replace==T) message("Duplicate genes will be averaged and merged to keep row.names unique")

    genes_from <- if (is.null(colGene)) rownames(dataIn) else dataIn[[colGene]]
    idx_match <- match(genes_from, df_mapping[[fromMapCol]])
    genes_to <- df_mapping[[toMapCol]][idx_match]
  
    if (replace) { # remove NAs
      if (is.null(colGene)) {
        # average identical gene names to ensure unique row names 
        dataIn_aggr <- aggregate(dataIn, by= list(genes_to), FUN=mean, na.rm=T)
        rownames(dataIn_aggr) <- dataIn_aggr[["Group.1"]]
        dataIn <- within(dataIn_aggr, rm("Group.1"))
      } else {
        if (na.rm) {
          dataIn <- dataIn[!is.na(genes_to),]
          genes_to <- genes_to[!is.na(genes_to)]
        }
        dataIn[[colGene]] <- genes_to
        colnames(dataIn)[which(colnames(dataIn)==colGene)] <- to
      }
    }  else {
      if (na.rm) {
        dataIn <- dataIn[!is.na(genes_to),,drop=F]
        dataIn[[to]] <- genes_to[!is.na(genes_to)]
      } else {
      dataIn[[to]] <- genes_to  
      }
    }
  } else if (class(dataIn)=="list") {
    dataIn <- lapply(dataIn, function(eachVec) {
      oldNames <- if (class(eachVec)== "numeric") {names(eachVec)} else if (class(eachVec)=="character") {eachVec}
      newNames <- df_mapping[[toMapCol]][match(oldNames, df_mapping[[fromMapCol]])]
      if (na.rm) {newNames <- newNames[!is.na(newNames)]
      eachVec <- eachVec[!is.na(newNames)]}
      if (class(eachVec)=="numeric") names(eachVec) <- newNames else eachVec <- newNames
      return(eachVec)
    })
  } else if (class(dataIn)=="numeric") {
    newNames <- df_mapping[[toMapCol]][match(names(dataIn), df_mapping[[fromMapCol]])]
    if (na.rm) {
      newNames <- newNames[!is.na(newNames)]
      dataIn <- dataIn[!is.na(newNames)]
      }
    names(dataIn) <- newNames
  } else if (class(dataIn)=="character") {
    dataIn <- df_mapping[[toMapCol]][match(dataIn, df_mapping[[fromMapCol]])]
    if (na.rm) dataIn <- dataIn[!is.na(dataIn)]
  }
  return(dataIn)
}

######################################################################
############################ GSEA homemade ###########################
######################################################################

GSEA_homemade <- function(datExpr,
                          C,
                          S,
                          p=1) {
  #' @usage compute the enrichment score using the GSEA algorithm 
  #' See Subramanian and Tamayo, 2005, PNAS
  #' @param datExpr = # N genes * k sample expression matrix 
  #' @param C metadata vector of length N with gene names
  #' @param S geneset: unordered vector of gene names associated with a pathway 
  #' @param p exponent to control how much gene correlation weights matter: 0 for equal weight, 1 for no change
  
  stopifnot(p %in% c(0,1))
  
  N = nrow(E)
  k = ncol(E)
  
  # correlate gene profiles with metadata
  R = abs(cor(t(datExpr), C)) %>% sort(.,decreasing=T) 
  # Evaluate the fraction of genes in S (‘‘hits’’) weighted by their correlation up to
  # position (gene) i 
  P_hit_S_i = function(i) sum(R[1:i][names(R)[1:i] %in% S]^p) / sum(R[names(R)[1:i] %in% S]^p)
  # Evaluate the fraction of genes not in S (‘‘misses’’) up to position (gene) i 
  P_miss_S_i = function(i) 1 / (N - sum((names(R) %in% S)[1:i])) 
  # Find the max difference in this running score
  runningscore = sapply(1:N, function(i) P_hit_S_i - P_miss_S_i)
  ES = runningscore[which.max(runningscore %>% abs)]
  
  # Compute significance by permuting C
  # TODO 
  
  # Find leading edge subset of geneset
  return(ES)
}

######################################################################
########################## TIMSHEL FUNCTIONS #########################
######################################################################

### Updated to Seurat 3.0, Jonatan Thompson 20190315

# from /projects/timshel/git/perslab-sc-library/seurat_functions/utils_seurat.R
filter_genes_in_seurat_data.dataframe <- function(df, 
                                                  seurat_obj, 
                                                  assay="RNA",
                                                  slot="data",
                                                  colname_gene, 
                                                  do.print=T) {
  #' @usage        filter data frame to contain only genes in seurat data.
  ### INPUT
  # ....
  ### OUTPUT
  # df            a data frame, filtered to contain only genes in the seurat data. All columns are returned.
  df <- as.data.frame(df) # make sure we start off with a data frame (and not a tibble).
  genes_seurat_data <- rownames(seurat_obj) # vector of genes
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
  genes_seurat_data <- rownames(seurat_obj) # vector of genes
  bool.genes_in_seurat_data <- x %in% genes_seurat_data # boolean
  if (do.print) {
    genes_not_found <- x[!bool.genes_in_seurat_data]
    print(sprintf("Filtered vector. Number of genes not found in Seurat data: %s", length(genes_not_found)))
    print(genes_not_found)
  }
  x.filtered <- x[bool.genes_in_seurat_data] # subset
  return(x.filtered)
}

# dotplot_timshel
# Pascal Timshel EDIT of Seurat code (3.0)
# Jonatan Thompson EDIT of Timshel edit (to Seurat 3.0)

######################## DOCUMENTATION #############################
# df.marker.panel.to.plot MUST contain columns "cell_type" and "gene_name"


################### MAKING FUNCTIONS AVAILABLE ###################
PercentAbove <- Seurat:::PercentAbove
# REF1: https://stackoverflow.com/questions/12178830/change-internal-function-of-a-package
# REF2: https://stackoverflow.com/questions/8743390/how-do-i-override-a-non-visible-function-in-the-package-namespace

### SEE ALSO:
# ?assignInNamespace
# ?fixInNamespace
# path.package("Seurat")

### data.to.plot snippet (seurat makes this)
# id	genes.plot	avg.exp	pct.exp	avg.exp.scale	
# 1	0	Acta2	7.31E-03	0.002913564	-0.208472888
# 2	0	Agrp	2.94E-02	0.011330528	-0.31387005
# 3	0	Aif1	8.09E-02	0.033667854	-0.283582648
# 4	0	Anxa2	6.67E-02	0.032696666	-0.801313835
# 5	0	Ccdc153	6.03E-02	0.021689867	-0.288510969
# 6	0	Ccl7	8.78E-04	0.000323729	-0.308715642

#################### START SCRIPT ########################

#' Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level of
#' 'expressing' cells (blue is high).
#'
#' @param object Seurat object
#' @param assay Seurat object assay
#' @param do.scale display mean expression scaled around mean by sd? If F (default), plot log1p(mean.expr)
#' @param genes.plot Input vector of genes
#' @param cols.use colors to plot
#' @param col.min Minimum scaled average expression threshold (everything smaller
#'  will be set to this)
#' @param col.max Maximum scaled average expression threshold (everything larger
#' will be set to this)
#' @param dot.min The fraction of cells at which to draw the smallest dot
#' (default is 0.05). All cell groups with less than this expressing the given
#' gene will have no dot drawn.
#' @param dot.scale Scale the size of the points, similar to cex
#' @param group.by Factor to group the cells by
#' @param plot.legend plots the legends
#' @param x.lab.rot Rotate x-axis labels
#' @param do.return Return ggplot2 object
#'
#' @return default, no return, only graphical output. If do.return=TRUE, returns a ggplot2 object
#'
#' @importFrom tidyr gather
#' @importFrom dplyr %>% group_by summarize_each mutate ungroup
#'
#' @export
#'
#' @examples
#' cd_genes <- c("CD247", "CD3E", "CD9")
#' DotPlot(object = pbmc_small, genes.plot = cd_genes)
#'

DotPlot_timshel <- function(
  object,
  assay="RNA",
  do.scale=F,
  genes.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  group.order=NULL,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = T,
  coord_flip = F,
  df.marker.panel.to.plot=NULL
  
) {
  
  
  # PT added: filter genes, to avoid errors with genes missing
  if (!is.null(df.marker.panel.to.plot)) {
    df.marker.panel.to.plot <- filter_genes_in_seurat_data.dataframe(df.marker.panel.to.plot, 
                                                                      seurat_obj=object, 
                                                                      colname_gene="gene_name", 
                                                                      do.print=T)
  }
  genes.plot <- filter_genes_in_seurat_data.vector(genes.plot, seurat_obj=object, do.print=F)
  
  if (! missing(x = group.by)) {
    Idents(object) <- object[[group.by]]
  }
  
  # Use the data slot because some cells may have been filtered out which are still present in the counts slot
  DefaultAssay(object) <- assay
  data.to.plot <- FetchData(object = object, slot = "data", vars=genes.plot)
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- Idents(object)
  
  data.to.plot %>% gather( 
    key = genes.plot,
    value = expression,
    -c(cell, id) # so 'gather' the genes.plot columns into a single column called expression. Don't touch cell and id columns
  ) -> data.to.plot # tidyr::gather takes multiple columns and collapses into key-value pairs, 
  # duplicating all other columns as needed.
  
  data.to.plot %>%
    group_by(id, genes.plot) %>% 
    #Most data operations are done on groups defined by variables. 
    # group_by() takes an existing tbl and converts it into a grouped tbl where operations are performed "by group". 
    # the output of group_by is a tbl where the variables that were grouped - in this case cell barcodes, which were
    # grouped by celltype ID - are missing 
    # ungroup() removes grouping
    # summarise replaces expression 
    summarize(
      avg.exp = mean(expm1(x = expression)), # convert back to raw counts in order to take mean
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  
  if (do.scale) {
    data.to.plot %>%
      ungroup() %>%
      group_by(genes.plot) %>%
      mutate(avg.exp.scale = scale(x = avg.exp)) %>%
      mutate(avg.exp.scale = MinMax(
        data = avg.exp.scale,
        max = col.max,
        min = col.min #This seems completely redundant..
      )) ->  data.to.plot
    
    data.to.plot$avg.exp <- data.to.plot$avg.exp.scale
  } else {
    # convert mean expression to log
    data.to.plot$avg.exp <- log1p(data.to.plot$avg.exp)
  }
  
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  
  data.to.plot$id <- factor(data.to.plot$id, levels=rev(group.order)) #reorder(x=data.to.plot$id, X=match(data.to.plot$id, group.order), FU)
 
  ### PT
  if (!is.null(df.marker.panel.to.plot)) {
    df.marker.panel.to.plot <- df.marker.panel.to.plot %>% 
      group_by(gene_name) %>% mutate(id = 0:(n()-1)) # OBS: column name "id" is also used in data.to.plot
    # ^ Create a sequential number (counter) for rows within each group of a dataframe 
    # ^ REF: https://stackoverflow.com/questions/11996135/create-a-sequential-number-counter-for-rows-within-each-group-of-a-dataframe
    df.marker.panel.to.plot <- df.marker.panel.to.plot %>% rename(genes.plot=gene_name)
    # order levels and clean gene names
    df.marker.panel.to.plot$genes.plot <- factor(
      x = sub(pattern = "-", replacement = ".", x = df.marker.panel.to.plot$genes.plot),
      levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
    )
  }
  
  ###
  p_dot <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradientn(colours=cols.use) +
    #scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) # <-- ORIGINAL
  
  if (! plot.legend) {
    p_dot <- p_dot + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p_dot <- p_dot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  
  if (coord_flip) p_dot <- p_dot + coord_flip()
  
  if (!is.null(df.marker.panel.to.plot)) {
    ### Marker plot
    p_markers <- ggplot(df.marker.panel.to.plot, aes(x=genes.plot)) + 
      geom_bar(aes(fill=cell_type), show.legend=F, color="white") +
      geom_text(aes(y=id, label=cell_type), angle=90, size=rel(3)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    # hjust=0 --> semi does not work
    # hjust="outward" --> does not work
    # hjust="inward" --> does not work
    # Inward always aligns text towards the center, and outward aligns it away from the center

    # Using cowplots package to align plots and axis 
    # REF: http://htmlpreview.github.io/?https://github.com/wilkelab/cowplot/blob/master/inst/doc/introduction.html
    p_dot <- p_dot + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    p_markers <- p_markers + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    p_grid <- plot_grid(p_dot, p_markers, ncol=1, align = "v", axis="lr", rel_heights=c(3,3))
    print(p_grid)
    # suppressWarnings(print(p))
  }
  
  print(p_dot)
  
  if (do.return) {
    return(p_dot)
    if (!is.null(df.marker.panel.to.plot)){
      return(p_grid)
      return(df.marker.panel.to.plot)
    }
  } else {
    print("Returning nothing")
    }
}