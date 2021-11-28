
#' Plot UMAP with specific styling
#'
#' @importFrom scater ggcells
#' @importFrom ggplot2 aes_string
#' @importFrom ggrastr geom_point_rast
#'
#' @param sce SingleCellExperiment object
#' @param colour_by The colData observation used to colour points.
#' @param point_size Size of the scatter points.
#' @param legend_position Position of the colour_by legend. Default is "right".
#' @param celltype_colours List of colours to map data values to. The values
#' will be matched in order.
#' @export
plotUMAP <- function(sce,colour_by="celltype",point_size=1,legend_position="right",
                     colours,alpha=1) {
  colnames(reducedDim(sce,"UMAP")) <- c("UMAP.1","UMAP.2")
  scater::ggcells(sce, aes_string(x="UMAP.1",y="UMAP.2",
                                  colour=as.factor(sce[[colour_by]]))) +
    scale_color_manual(values = colours, name = "") +
    xlab("UMAP1") + ylab("UMAP2") +
    ggrastr::geom_point_rast(size = point_size,alpha=alpha) + theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position=legend_position,
          aspect.ratio = 1) +
    guides(col = guide_legend(override.aes = list(size = 5))) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.grid = element_blank())

}


#' Plot dimensionality reduction highlighting observations
#'
#' @param sce SingleCellExperiment object
#' @param dimred A string indicating the reduced dimesnion result in
#' `reducedDims(sce)` to plot.
#' @param colour_by A string indicating the `colData(sce)` column to
#' colour points by.
#' @param highlight A list of strings specifying the `colour_by` values to
#' highlight.
#' @param colours A list of colours of colour palette passed to `scale_color_manual`.
#' @param legend_position String passed to legend_position paramater of ggplot theme.
#' @param legend_title Boolean indicating whether to plot a legend title.
#'
#' @export
plotDimredHighlight <- function(sce,dimred="TSNE",colour_by="AO_annotation",
                              highlight=c("Serotonergic neurons",
                                          "Early AP", "Late AP",
                                          "Other neural"),
                              colours=echino.AO_colours,point_size=1,
                              plot_order=NULL,
                              legend_position="right",
                              legend_title=T) {
  df <- reducedDim(sce,dimred)
  df <- as.data.frame(df)
  colnames(df) <- c(paste0(dimred,"_1"),paste0(dimred,"_2"))
  df <- cbind(df,sce[[colour_by]])
  colnames(df)[3] <- colour_by

  if(!is.null(plot_order)) {
    df <- df[plot_order,]
  }

  df_highlight <- df[sce[[colour_by]]%in%highlight,]

  # order of plot according to order of highlight
  df_highlight[,colour_by] <- factor(df_highlight[,colour_by],levels=highlight)
  df_highlight <- df_highlight[order(df_highlight[,colour_by],decreasing = T),]

  p <- ggplot(df,aes_string(x=paste0(dimred,"_1"),y=paste0(dimred,"_2"))) +
    ggrastr::geom_point_rast(data=df,alpha=0.8,size=point_size,shape=16,color="grey")
  p <- p+
    ggrastr::geom_point_rast(data=df_highlight,
               aes_string(x=paste0(dimred,"_1"),y=paste0(dimred,"_2"),
                          color=colour_by),size=point_size,shape=16)

  if(!is.null(colours)) {
    p <- p + scale_color_manual(values = colours, name = "")
  }

  legend_title <- ifelse(legend_title,colour_by,"")

  p <- p + theme_classic() +

    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = legend_position) +
    guides(col = guide_legend(override.aes = list(size = 5),
                              title=legend_title)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.grid = element_blank())


  return(p)
}





getGeneMeansFast <- function(gene, mat, celltypes){
  expr <- mat[gene,]
  medians = aggregate(formula = expr ~ celltypes, FUN = mean)
  out = matrix(medians$expr, nrow = 1, dimnames = list(gene, medians$celltypes))
  return(out)
}

getExpressionFast <- function(sce,genes,group_by="celltype",norm="min-max") {
  genes <- genes[!is.na(genes)]
  mat <- logcounts(sce)[genes,]
  means <- lapply(genes, getGeneMeansFast, mat = mat, celltypes = sce[[group_by]])
  combined <- do.call(rbind,means)
  celltypes <- colnames(combined)
  if(norm=="min-max") {
    combined <- sweep(combined, 1, apply(combined, 1, max), "/")
  } else if(norm=="z-scale") {
    combined <- t(apply(combined, 1, scale))
    colnames(combined) <- celltypes
  }
  return(combined)
}
