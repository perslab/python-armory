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
  if (!is.null(genes)) if (any(!genes %in% rownames(seurat_obj))) {
    stop(paste0(paste0(genes[!genes %in% rownames(seurat_obj)], collapse= " "), " not found in seurat object raw data"))
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
    list_idx_gene <- lapply(X=comb, FUN=function(gene) GetAssayData(seurat_obj, slot="data")[gene, ]>= detection.threshold)
    idx_comb_intersect <- Reduce(x=list_idx_gene, f = '&')
    
    # Get the **union** of the idx of the cells expressing the genes which are **not** in the combination (the complementary set of genes)
    if (length(comb_complement)>= detection.threshold) {
      list_idx_gene_complement <- lapply(X=comb_complement, FUN=function(gene) GetAssayData(seurat_obj, slot="data")[gene, ]>= detection.threshold)
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
  vec_labels <- character(length(Idents(seurat_obj)))
  
  for (i in 1:length(list_idx_comb)) {
    vec_labels[list_idx_comb[[i]]] <- names_comb[i]
  }
  
  vec_labels[nchar(vec_labels)==0] <- "None"#NA_character_
  
  # Add vector of labels to Seurat object as metadata
  df_labels <- data.frame(feats_to_plot = vec_labels, row.names = colnames(seurat_obj))
  
  seurat_obj <- AddMetaData(object=seurat_obj, metadata=df_labels)
  
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
  list_args[["label"]] = F
  list_args[["cols"]] = colors_plot_hex
  list_args[["do.return"]] = T
  list_args[["group.by"]] = "feats_to_plot"
  list_args[["no.legend"]] = F
  list_args[["plot.title"]] = plot.title
  list_args[["na.value"]] = "grey90"
  
  # Make the plot
  p <- do.call(what = DimPlot, 
               args= list_args)
  
  return(p)
}

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
  # ortholog mapping:
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
        vec_logicalNA <- is.na(genes_to)
        dataIn <- dataIn[!vec_logicalNA,,drop=F]
        dataIn[[to]] <- genes_to[!vec_logicalNA]
      } else {
      dataIn[[to]] <- genes_to  
      }
    }
  } else if (class(dataIn)=="list") {
    dataIn <- lapply(dataIn, function(eachVec) {
      oldNames <- if (class(eachVec)== "numeric") {
        names(eachVec)
      } else if (class(eachVec)=="character") {
          eachVec
        }
      newNames <- df_mapping[[toMapCol]][match(oldNames, df_mapping[[fromMapCol]])]
      if (na.rm) {
        vec_logicalNA <- is.na(newNames)
        eachVec <- eachVec[!vec_logicalNA]
        newNames <- newNames[!vec_logicalNA]
      }
      if (class(eachVec)=="numeric") names(eachVec) <- newNames else eachVec <- newNames
      return(eachVec)
    })
  } else if (class(dataIn)=="numeric") {
    newNames <- df_mapping[[toMapCol]][match(names(dataIn), df_mapping[[fromMapCol]])]
    if (na.rm) {
      vec_logicalNA <- is.na(newNames)
      newNames <- newNames[!vec_logicalNA]
      dataIn <- dataIn[!vec_logicalNA]
      }
    names(dataIn) <- newNames
  } else if (class(dataIn)=="character") {
    dataIn <- df_mapping[[toMapCol]][match(dataIn, df_mapping[[fromMapCol]])]
    if (na.rm) dataIn <- dataIn[!is.na(dataIn)]
  }
  return(dataIn)
}

######################################################################
############################ GSEA homemade workflow ###########################
######################################################################



if (F) {
  
  # Generate null replicas of the gene weight vector (assigning a value to all genes) used in GSEA
  require(boot)
  
  path_data = ""
  data = load_obj(path_data) # add dataset here
  statistic = mystatisticfunction # A function which when applied to data returns a vector 
  # containing the statistic(s) of interest. When sim = "parametric", the first argument 
  # to statistic must be the data. For each replicate a simulated dataset returned by ran.gen 
  # will be passed. In all other cases statistic must take at least two arguments. 
  # The first argument passed will always be the original data. The second will be a vector of indices, 
  # frequencies or weights which define the bootstrap sample.
  R = 10000 # number of bootstrap reps
  indices = myindices # a vector of 'case' indices
  parallelOperation = "multicore"
  options("boot.parallel"="multicore")
  
  randomSeed = 12345
  sim = "permutation" # i.e. draw full sample without replacement
  
  bootout <- boot::boot(data=data,
                        statistic=statistic,
                        R = R,
                        seed=randomSeed,
                        sim=sim,
                        stype="i",
                        parallel=parallelOperation,
                        ncpus=30,
                        cl=NULL)
  
}


fnc_GSEAperslab <- function(list_vec_S,
                          fnc_geneScore=NULL,
                          vec_R = NULL,
                          list_vec_Rnull = NULL,
                          vec_geneNames = NULL,
                          negEnrichment=F,
                          p=1,
                          nRep=1000,
                          randomSeed=12345,
                          timeout = 43200,
                          ...)  {
  #' See Subramanian and Tamayo, 2005, PNAS
  #' @usage compute the enrichment score using the GSEA algorithm. 
  #' @param list_vec_S named list of unordered geneset vectors; gene naming should match datExpr 
  #' @param fnc_geneScore function to score genes, e.g. correlation or diff expr. Class function. default NULL
  #' @param vec_R precomputed named gene score vector. Default NULL
  #' @param list_vec_Rnull list of precomputed NULL named gene score vectors. Default NULL
  #' @param vec_geneNames vector of geneNames in the order they are given 
  #' @param negEnrichment take into account negative enrichment (depletion) when computing p-values (i.e. two tailed)
  #' @param p exponent to control how much gene correlation weights matter: 0 for equal weight, 1 for no change, defaults to 1
  #' @param nRep number of random permutations, defaults to 1000
  #' @param randomSeed random seed
  #' @param timeout number of seconds to allow for computing ES, defaults to 43200 (12 hours)
  #' @param ... appropriately named and formatted arguments to pass to fnc_geneScore. 
  #' The first must be a numeric or logical vector of binary sample annotations indicating condition vs control
  #' @return a list containing 
  #'           vec_E: a named vector of E scores 
  #'           vec_p.value: a named vector or unadjusted p-values
  #'           parameter values of 
  #'               fnc_geneScore, 
  #'               negEnrichment, 
  #'               p, 
  #'               nRep, 
  #'               randomSeed
  #'           vec_R: gene score vector (either as provided or computed)
  #'           list_vec_Rnull: list of named NULL gene score vectors (either as provided or computed)
  #'           mat_ESnull: matrix of null ES scores. dimensions nRep * length(list_vec_S)
  #' @depends 
  #' parallel package
  #' utility function safeParallel

  require(parallel)
  require(magrittr)
  
  set.seed(seed = randomSeed)
  
  # check inputs
  stopifnot(p %in% c(0,1), 
            !is.null(names(list_vec_S))) # need to have two levels
  if (!is.null(list_vec_Rnull)) if(length(list_vec_Rnull)<1000) warning("list_vec_Rnull contains fewer than 1000 null replicates")
  stopifnot(!is.null(fnc_geneScore) | !is.null(list_vec_Rnull) & !is.null(vec_R))
  if (!is.null(fnc_geneScore) & is.null(vec_geneNames)) stop("fnc_geneScore needs the argument vec_geneNames")
  
  # take ... args
  list_args <- list(...)

  # If required convert fnc_geneScore from object to character
  # if (!"character" %in% class(fnc_geneScore)) {
  #   fnc_geneScore <- as.character(quote(fnc_geneScore))
  # }

  # compute scores for all genes in expression data using the score function
  if (is.null(vec_R)) {
    fnc_geneScoreWrap <- function(fnc_geneScore, list_args)
      {
      vec_R <- do.call(what=fnc_geneScore, args=list_args) 
      if(is.null(names(vec_R))) names(vec_R) <- vec_geneNames
    }
    # Compute named vector of gene scores
    vec_R <- fnc_geneScoreWrap(fnc_geneScore=fnc_geneScore, list_args=list_args)
  } 
  
  vec_R <- sort(vec_R, decreasing=T) 

  # compute ES scores
  # see https://www.pnas.org/content/pnas/102/43/15545.full.pdf
  fnc_ES <- function(vec_R, vec_S, p) {
    
    vec_R <- sort(vec_R, decreasing=T)
    # First check if there is any overlap at all
    if (sum(names(vec_R) %in% vec_S)==0) return(0)
    
    vec_logicalRinS <- names(vec_R) %in% vec_S
    # the cum sum over the elements of vec_R which are in vec_S, divided by the sum
    vec_Rhit <- vec_R
    vec_Rhit[!vec_logicalRinS] <- 0
    vec_hitCumSum <- vec_Rhit %>% abs %>% '^'(p) %>% cumsum
    vec_Phit <- vec_hitCumSum/vec_hitCumSum[length(vec_hitCumSum)]
    # nb: the denominator is just the sum
    # cum sum over elements of vec_R not in S, divided by the sum
    vec_missCumSum <- vec_logicalRinS %>% '!'(.) %>% as.numeric %>% cumsum 
    vec_Pmiss = vec_missCumSum / vec_missCumSum[length(vec_missCumSum)]
    # nb: the denominator is just the sum
    
    # return max divergence
    max(vec_Phit-vec_Pmiss)
  }
  
  vec_ES <- sapply(list_vec_S, function(vec_S) fnc_ES(vec_S=vec_S,
                                              vec_R=vec_R, 
                                              p=p), simplify=T) 
  
  # get rid of logical(0)s
  vec_ES[is.na(vec_ES==0)] <- 0
  
  # Compute significance by permuting labels
  ## Make list of nRep random gene score vectors using permuted labels C
  
  if (is.null(list_vec_Rnull)) {
    message("Preparing NULL gene weight vectors")
    list_vec_Rnull <- lapply(1:nRep, function(i){
      list_argsCRand <- list_args
      list_argsCRand[[1]] <- sample(x=list_args[[1]], size = length(list_args[[1]]), replace = F)
      suppressMessages({fnc_geneScoreWrap(fnc_geneScore=fnc_geneScore,
                        list_args = list_argsCRand)})
    })
  } else {
    # if NULL gene score vectors have been precomputed, sort them
    list_vec_Rnull <- lapply(list_vec_Rnull, sort, decreasing=TRUE)
    nRep = length(list_vec_Rnull)
  }

  # Compute ES scores for the 'random' R gene weight vectors
  message("Computing GSEA E scores for null gene weight vectors")
  mat_ESnull <- sapply(list_vec_S, function(vec_S) {
    suppressMessages({
      safeParallel(fun=fnc_ES, 
                   list_iterable=list("vec_R"=list_vec_Rnull),
                   timeout = timeout, 
                   vec_S = vec_S,
                   p = p,
                   simplify = T)})
    }, simplify=T)
  
  mat_ESnull[is.na(mat_ESnull==0)] <- 0
  
  # Compute empirical statistical significance
  message("Computing empirical p-values")
  
  fnc_pval <- function(i, ES){
    vec_ESnull <- mat_ESnull[,i]
    if (negEnrichment) sum(abs(c(vec_ESnull,ES))>=abs(ES)) else sum(c(vec_ESnull,ES)>=ES) %>% 
      '/'(length(c(vec_ESnull,ES)))
      # count at both tails using absolute values
      # note that multiplying the p-value for the one-tailed test would be assuming 
      # that the null distribution is symmetrical, which it probably isn't
      # see https://www.omicsonline.org/open-access/twotailed-pvalues-calculation-in-permutationbased-tests-a-warning-against-asymptotic-bias-in-randomized-clinical-trials-2167-0870.1000145.php?aid=19520
      # see https://stats.stackexchange.com/questions/140107/p-value-in-a-two-tail-test-with-asymmetric-null-distribution
  }
  
  list_iterable <- list("i"=1:length(list_vec_S), 
                        "ES"=vec_ES)
    
  vec_P <- safeParallel(fun = fnc_pval,
                        list_iterable=list_iterable,
                        SIMPLIFY=T)
  
  names(vec_P) <- names(list_vec_S)
  
  # prepare outputs
  list_out <- list("vec_ES" = vec_ES, 
                  "vec_p.value" = vec_P,
                  "fnc_geneScore" = if (!is.null(fnc_geneScore)) as.character(quote(fnc_geneScore)) else NULL,
                  "negEnrichment"= negEnrichment,
                  "p" = p,
                  "nRep" = nRep,
                  "randomSeed" = randomSeed,
                  "vec_R" = vec_R,
                  "list_vec_Rnull"=list_vec_Rnull,
                  "mat_ESnull"=mat_ESnull)
                 
  return(list_out)
}


# tests

if (F) {
  
  
  ### TEST 1 ###
  # No signal expected

  list_genesets <- load_obj("/projects/jonatan/genesets/c2.all.mmusculus.symbols.RData")
  list_vec_S <- list_genesets[1:10] %>% lapply(X=., FUN=function(vec_S) {
    paste0(substr(vec_S,1,1), tolower(substr(x=vec_S, 2, length(vec_S))))
  })

  seuratChen <- readRDS("/projects/jonatan/data/Chen_all_cells_seu3.RDS.gz")
  seuratChen <- FindVariableFeatures(seuratChen, nfeatures = 2000)
  seuratChenSub <- subset(seuratChen, subset=cell_type %in% c("Ependy", "GABA1"), 
                          features=VariableFeatures(seuratChen))

  
  list_args=list()
  list_args[["x"]] <- GetAssayData(object=seuratChenSub, slot="data") %>% as.matrix %>% t
  list_args[["y"]] <- seuratChenSub$cell_type %>% as.factor %>% as.numeric %>% '-'(1) %>% as.matrix(.,nrow=length(.))
  rownames(list_args[["y"]]) <- seuratChenSub$cell_type
  list_args[["method"]]= "pearson"

  fnc_geneScore="cor"
  negEnrichment=F
  p=1
  nRep=10
  randomSeed=12345
  
  ### TEST 2 ###
  # expression data and annotations
  mat_counts <- load_obj("/raid5/home/cbmr/lhv464/astrogsea.RDS")
  mat_counts <- mat_counts[,grepl("FGF1|PF", colnames(mat_counts))]
  keep <- rowSums(mat_counts >= 10) > 10
  mat_counts <- mat_counts[keep,]
  trt<-as.factor(sapply(strsplit(sapply(strsplit(colnames(mat_counts),"_"),"[",2),"\\."),"[",1))
  
  # genesets
  list_genesets <- load_obj("/projects/jonatan/genesets/bmi_brain_hsapiens_190408.RDS")
  df_orthologMapping <- load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz")
  df_ensemblMapping <- load_obj("/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz")
  list_vec_S <- sample(list_genesets, size=5, replace=F)
  list_vec_S <- lapply(list_vec_S, function(geneset) {
    sapply(geneset, function(gene) {
      paste0(toupper(substr(gene,1,1)), tolower(substr(gene,2,nchar(gene))))
      })
    })
  
  list_args <- list("trt"=trt, "mat_counts"=mat_counts)
  
  fnc_DESeq2Wrapper <- function(trt, mat_counts) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat_counts,
                                          colData = data.frame("trt"=trt),
                                          design = ~ 0 + trt)  
    dds<-DESeq2::DESeq(dds)
    resultsObj <- DESeq2::results(dds, contrast=c("trt", "FGF1", "PF"))
    return(resultsObj$log2FoldChange)
  }
  fnc_geneScore = fnc_DESeq2Wrapper
  vec_geneNames = rownames(mat_counts)
  negEnrichment=F
  p=1
  nRep=10
  randomSeed=12345
  

}



fnc_GSEABroad <- function(gsea_jar="/projects/jonatan/tools/gene_set_enrichment/GSEA/gsea-3.0.jar",
                        dest_gmt_file,
                        rnk_file,
                        num_randomizations,
                        analysis_name,
                        rand_working_dir) {
  #' @usage call GSEA java implementation
  #' @param gsea_jar 
  #' @param dest_gmt_file
  #' @param rkn_file
  #' @param num_randomizations a 
  #' @param analysis_name
  #' @param rand_working_dir
  #' @return NULL
  
  # copied from https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/Supplemental_protocol_4_manual_phenotype_rand_with_edgeR.html#118_run_gsea
  
  start_gs_perm <- Sys.time()
  command <- paste("java  -Xmx1G -cp ",gsea_jar,  " xtools.gsea.GseaPreranked -gmx ", 
                   dest_gmt_file, "-rnk " ,rnk_file, 
                   "-collapse false -nperm ",num_randomizations, 
                   " -permute gene_set -scoring_scheme weighted -rpt_label ",
                   paste(analysis_name,"gsrand",sep="_"),
                   " -num 100 -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out " ,
                   rand_working_dir, "-gui false > gsea_output.txt",sep=" ")
  system(command)
  stop_gs_perm <- Sys.time()
  difftime(stop_gs_perm,start_gs_perm,"mins")
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
