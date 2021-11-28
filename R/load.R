
library(Matrix)
library(SingleCellExperiment)
library(scater)
library(scran)
library(data.table)


foster.path <- "G:/My Drive/Postgrad/PhD/Projects/data/Foster_2020/"

loadFoster2020 <- function(foster.path) {

  # Load SingleCellExperiment
  echino.sce <- readRDS(paste0(foster.path,"sce.rds"))

  # Reannotate stage
  echino.sce$stage <- ""
  echino.sce$stage[echino.sce$orig.ident=="Sp1"] <- "8-cell"
  echino.sce$stage[echino.sce$orig.ident=="Sp2"] <- "64-cell"
  echino.sce$stage[echino.sce$orig.ident=="Sp3"] <- "Morula"
  echino.sce$stage[echino.sce$orig.ident=="EB"] <- "Early blastula"
  echino.sce$stage[echino.sce$orig.ident=="HB"] <- "Hatched blastula"
  echino.sce$stage[echino.sce$orig.ident=="MB"] <- "Mesenchyme blastula"
  echino.sce$stage[echino.sce$orig.ident=="EG"] <- "Early gastrula"
  echino.sce$stage[echino.sce$orig.ident=="LG"] <- "Late gastrula"
  echino.sce$stage <- factor(echino.sce$stage,
                             levels=c("8-cell","64-cell","Morula","Early blastula",
                                      "Hatched blastula","Mesenchyme blastula",
                                      "Early gastrula","Late gastrula"))

  # Label cell types from Seurat clusters
  echino.sce$original_annotation <- "Other"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(3,9,17,18)] <- "Neural"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(20)] <- "Germline"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(1,4,12,13)] <- "Oral ectoderm"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(16,19)] <- "Skeleton"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(6,8,14)] <- "Endoderm"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(2,7)] <- "Ciliated cells"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(11)] <- "Pigment cells"
  echino.sce$original_annotation[echino.sce$seurat_clusters %in% c(0,5)] <- "Aboral ectoderm"

  # Add refined AO annotations
  AO_annotations <- read.csv(paste0(foster.path,"AO_annotations_v2.tsv"),
                                  sep="\t",header = T, row.names=1)
  echino.sce$AO_annotation <- AO_annotations[,1]


  # Load UMAP
  umap <- read.table(paste0(foster.path,"umap.tsv"),sep="\t", header=T)
  rownames(umap) <- colnames(echino.sce)
  reducedDim(echino.sce, "UMAP") <- umap

  # Load PCA
  pca <- data.table::fread(paste0(foster.path,"echino_fastMNN_pcs.tsv"),sep="\t",
                           header=F, data.table = F)[,-1]
  colnames(pca) <- paste0("PC_", 1:ncol(pca))
  rownames(pca) <- colnames(echino.sce)
  reducedDim(echino.sce,"PCA") <- pca


  # Edit annotations for heatmap plots
  echino.sce$AO_heatmap_annotation <- echino.sce$AO_annotation
  echino.sce$AO_heatmap_annotation[echino.sce$AO_heatmap_annotation=="Serotonergic neurons"] <- "SN"
  echino.sce$AO_heatmap_annotation[
    echino.sce$AO_heatmap_annotation%in%c("Oral ectoderm", "Aboral ectoderm",
                                          "Ciliated cells")
    ] <- "Other ectoderm"

  echino.sce$AO_heatmap_annotation[
    echino.sce$AO_heatmap_annotation%in%c("Skeleton", "Pigment cells")
  ] <- "Mesoderm"

  echino.sce$AO_heatmap_annotation[
    echino.sce$AO_heatmap_annotation%in%c("Germline", "Other")
  ] <- "Other"

  # Enforce specific cell type ordering
  echino.sce$AO_heatmap_annotation <- factor(echino.sce$AO_heatmap_annotation,
                                             levels=c("Early AP","Late AP",
                                                      "SN","Other neural",
                                                      "Other ectoderm",
                                                      "Mesoderm","Endoderm",
                                                      "Other"))

  echino.sce$AO_heatmap_stage <- echino.sce$stage
  levels(echino.sce$AO_heatmap_stage) <- c("8c","64c","M","EB","HB","MB","EG",
                                           "LG")

  return(echino.sce)
}




wagner.path <- "G:/My Drive/Postgrad/PhD/Projects/data/Wagner_2018/"


loadWagner2018 <- function(path) {

  meta <- data.frame(read.csv(file = paste0(path,'meta.csv'), header=TRUE))
  colnames(meta)[10:11] = c("celltype","stage")

  cell_ids = meta$unique_cell_id
  rownames(meta) = cell_ids
  genes <- data.frame(read.csv(file = paste0(path,'genes.csv'), header=TRUE))$index

  counts <- t(readMM(paste0(path,"matrix.mtx")))
  colnames(counts) = cell_ids
  rownames(counts) = genes

  corrected_pcs <- read.table(paste0(path,"corrected_pcs.tsv"),sep="\t")

  sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = meta,
                              reducedDims=list(PCA=corrected_pcs))

  size_factors <- read.csv(paste0(path,"size_factors.tsv"),row.names = 1)
  sce <- logNormCounts(sce,size_factors=size_factors$x)

  return(sce)

}


exportSample <- function(sce,ncells,out_path,out_name="sce_sample.RDS"){
  cells <- sample(seq(1,ncol(sce),by=1),ncells,replace=FALSE)
  sce_sample <- sce[,cells]
  saveRDS(sce_sample,paste0(out_path,out_name))
}



loadWagnerSpring <- function(wagner.path) {

  celltypes <- t(read.csv(paste0(wagner.path,"spring_data/diencephalon/cell_groupings.csv"),sep=",",header=F))
  colnames(celltypes) <- celltypes[1,]
  celltypes <- celltypes[2:nrow(celltypes),]
  rownames(celltypes) <- NULL

  cell_names <- read.csv(paste0(wagner.path,"spring_data/diencephalon/original_cell_indices.txt"),header = F,col.names=c("cell"))
  norm_counts <- read.csv(paste0(wagner.path,"spring_data/diencephalon/expr.csv.gz"),header = F,row.names=1)
  colnames(norm_counts) <- paste0("cell_",1:ncol(norm_counts))
  #coords <- read.csv(paste0(wagner.path,"spring_data/diencephalon/coordinates.csv"),header=F,row.names = 1)
  rownames(coords) <- paste0("cell_",1:nrow(coords))

  edges <- read.csv(paste0(wagner.path,"spring_data/"))

  sce <- SingleCellExperiment(assays=list(logcounts=as.matrix(norm_counts)),
                              colData=celltypes,
                              reducedDims = list(SPRING=coords))

  return(sce)
}


loadWagnerStich <- function() {
  norm_counts <- read.csv(paste0(wagner.path,"STITCH/export_directory/counts.csv"))
  edges <- read.csv(paste0(wagner.path,"STITCH/export_directory/edges.csv"),header=F,
                    sep=";")

}

reorderSTITCH <- function(zeb.sce,stitch_meta) {

  k <- 1
  zeb.sce$stitch_index <- -1

  for(i in c("24hpf","18hpf","14hpf","10hpf","8hpf","6hpf","4hpf")) {
    l <- k + sum(stitch_meta$NodeLabels==i)
    zeb.sce$stitch_index[zeb.sce$stage == i][stitch_meta$OriginalName[stitch_meta$NodeLabels==i]] = k:(l-1)
    k <- l
  }

  return(zeb.sce)

}




