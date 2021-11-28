
#' UMAP wrapper function
#'
#' Computes UMAP projection on principal components obtained from reducedDim slot.
#'
#' @param sce SingleCellExperiment object
#' @param min_dist The effective minimum distance between embedded points.
#' @param n_neighbours Integer scalar, number of nearest neighbors to identify
#'  when constructing the initial graph.
#' @param n_dims The dimension of the space to embed into.
#' @return SingleCellExperiment object with UMAP results in a reducedDim slot.
#'
#' @importFrom umap umap
#' @export
computeUMAP <- function(sce,min_dist=0.1,n_neighbours=50,n_dims=2) {
  out <- umap::umap(reducedDim(sce,"PCA"),
                    min_dist=min_dist,n_neighbors=n_neighbours,method="umap-learn",n_components=n_dims)
  reducedDim(sce,"UMAP") <- out$layout
  return(sce)
}




