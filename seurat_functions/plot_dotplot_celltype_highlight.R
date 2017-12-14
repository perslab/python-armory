# Pascal Timshel EDIT of Seurat code (2.0)
# October 2017.

library(Seurat)
library(cowplot)


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
  genes.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE,
  df.marker.panel.to.plot
) {
  
  
  # PT added: filter genes, to avoid errors with genes missing
  df.marker.panel.to.plot <- filter_genes_in_seurat_data.dataframe(df.marker.panel.to.plot, seurat_obj=object, colname_gene="gene_name", do.print=T)
  genes.plot <- filter_genes_in_seurat_data.vector(genes.plot, seurat_obj=object, do.print=F)
  
  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  data.to.plot$genes.plot <- factor(
    x = data.to.plot$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  
  ### PT
  df.marker.panel.to.plot <- df.marker.panel.to.plot %>% group_by(gene_name) %>% mutate(id = 0:(n()-1)) # OBS: column name "id" is also used in data.to.plot
  # ^ Create a sequential number (counter) for rows within each group of a dataframe 
  # ^ REF: https://stackoverflow.com/questions/11996135/create-a-sequential-number-counter-for-rows-within-each-group-of-a-dataframe
  df.marker.panel.to.plot <- df.marker.panel.to.plot %>% rename(genes.plot=gene_name)
  # order levels and clean gene names
  df.marker.panel.to.plot$genes.plot <- factor(
    x = sub(pattern = "-", replacement = ".", x = df.marker.panel.to.plot$genes.plot),
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )

  p_dot <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
    geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) # <-- ORIGINAL
  
  if (! plot.legend) {
    p_dot <- p_dot + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p_dot <- p_dot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  
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
  
  if (do.return) {
    # print("Returning nothing")
    return(p_grid)
    # return(df.marker.panel.to.plot)
  }
}






#################### START SCRIPT - FACET WRAP - NOT SO NICE ########################

#' Dot plot visualization
#'
#' Intuitive way of visualizing how gene expression changes across different
#' identity classes (clusters). The size of the dot encodes the percentage of
#' cells within a class, while the color encodes the AverageExpression level of
#' 'expressing' cells (blue is high).
#'
#' @param object Seurat object
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

DotPlot_timshel_facet_wrap <- function(
  object,
  genes.plot,
  cols.use = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by,
  plot.legend = FALSE,
  do.return = FALSE,
  x.lab.rot = FALSE,
  df.marker.panel.to.plot
) {
  if (! missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot %>% gather(
    key = genes.plot,
    value = expression,
    -c(cell, id)
  ) -> data.to.plot
  data.to.plot %>%
    group_by(id, genes.plot) %>%
    summarize(
      avg.exp = mean(expm1(x = expression)),
      pct.exp = PercentAbove(x = expression, threshold = 0)
    ) -> data.to.plot
  data.to.plot %>%
    ungroup() %>%
    group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>%
    mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = col.max,
      min = col.min
    )) ->  data.to.plot
  # data.to.plot$genes.plot <- factor(
  #   x = data.to.plot$genes.plot,
  #   levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  # ) # PT outcommented
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  
  ### PT
  df.marker.panel.to.plot <- df.marker.panel.to.plot %>% group_by(gene_name) %>% mutate(id = 0:(n()-1)) # OBS: column name "id" is also used in data.to.plot
  # ^ Create a sequential number (counter) for rows within each group of a dataframe 
  # ^ REF: https://stackoverflow.com/questions/11996135/create-a-sequential-number-counter-for-rows-within-each-group-of-a-dataframe
  df.marker.panel.to.plot <- df.marker.panel.to.plot %>% rename(genes.plot=gene_name)
  df.marker.panel.to.plot$panel <- "b"
  
  data.to.plot$panel <- "a"
  data.to.plot$id <- as.integer(data.to.plot$id) # TO AVOID ERROR -->  Error in bind_rows_(x, .id) :Column `id` can't be converted from integer to factor 
  data.combined <- bind_rows(df.marker.panel.to.plot, data.to.plot)
  
  # order levels and clean gene names
  data.combined$genes.plot <- factor(
    x = data.combined$genes.plot,
    levels = rev(x = sub(pattern = "-", replacement = ".", x = genes.plot))
  )
  
  # Align two plots on a page
  # REF: https://github.com/tidyverse/ggplot2/wiki/Align-two-plots-on-a-page
  p <- ggplot(data = data.combined, mapping = aes(x = genes.plot)) +
    facet_grid(panel~., scale="free") + 
    geom_point(data=data.to.plot, mapping = aes(y=as.factor(id), size = pct.exp, color = avg.exp.scale)) +
    geom_bar(data=df.marker.panel.to.plot, mapping = aes(fill=cell_type), show.legend=F, color="white") +
    geom_text(data=df.marker.panel.to.plot, mapping = aes(y=id, label=cell_type)) +
    scale_radius(range = c(0, dot.scale)) +
    scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  # p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, y = id)) +
  #   geom_point(mapping = aes(size = pct.exp, color = avg.exp.scale)) +
  #   scale_radius(range = c(0, dot.scale)) +
  #   scale_color_gradient(low = cols.use[1], high = cols.use[2]) +
  #   theme(axis.title.x = element_blank(), axis.title.y = element_blank()) # <-- ORIGINAL
  
  if (! plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(data.combined)
  }
}

