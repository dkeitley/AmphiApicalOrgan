
library(ApicalOrgan2021)

set.seed(44)


prepareEdgeR <- function(sce, group_by, group1, group2, block_by) {

  # Prepare groups
  if ((group1 != 'All') & (group2 != 'All')) {
    sce_sub <- sce[,colData(sce)[[group_by]] %in% c(group1, group2)]
    sce_sub$group <- colData(sce_sub)[[group_by]]
  } else if (group1 == 'All') {
    sce_sub <- sce
    sce_sub$group <- colData(sce_sub)[[group_by]]
    sce_sub$group[sce_sub$group != group2] <- 'All'
  } else if (group2 == 'All') {
    sce_sub <- sce
    sce_sub$group <- colData(sce_sub)[[group_by]]
    sce_sub$group[sce_sub$group != group1] <- 'All'
  }

  sce_sub$batch <- colData(sce_sub)[[block_by]]
  groups <- unique(sce_sub$group)
  sce_sub$group <- factor(sce_sub$group, levels = groups)

  return(sce_sub)
}



runEdgeR <- function(sce, min_detection_rate_per_group=0.1, min.logFC=2,
                     threshold_fdr=0.1) {


  # input should already have groups
  groups <- unique(sce$group)

  # calculate detection rate per gene
  cdr.dt <- data.table(
    gene = rownames(sce),
    detection_rate_A = rowMeans(counts(sce[,sce$group==groups[1]])>0),
    detection_rate_B = rowMeans(counts(sce[,sce$group==groups[2]])>0)
  ) %>% setnames(c("ens_id",sprintf("detection_rate_%s",groups[1]),sprintf("detection_rate_%s",groups[2])))


  # Filter genes by min detection rate per group
  cdr_A <- rowMeans(counts(sce[,sce$group==groups[1]])>0) >= min_detection_rate_per_group
  cdr_B <- rowMeans(counts(sce[,sce$group==groups[2]])>0) >= min_detection_rate_per_group
  sce <- sce[cdr_B | cdr_A,]

  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce, type="edgeR")

  # Define design matrix (with intercept)
  cdr <- colMeans(counts(sce)>0)
  design <- model.matrix(~ cdr + sce$batch + sce$group)

  # Estimate dispersions
  sce_edger  <- estimateDisp(sce_edger, design)

  # Fit GLM
  fit <- glmQLFit(sce_edger, design, coef=paste0('sce$group', groups[2]))

  # Test
  lrt <- glmQLFTest(fit, coef=paste0('sce$group', groups[2]))

  # Construct output data.frame
  out <- topTags(lrt, n=nrow(lrt))$table


  return(out)

}




###################
## Parse options ##
###################
# SCE path
# Group by variable
# Group 1
# Group 2
# Blocking factor
# Save path


option_list = list(
  make_option(c("-s", "--sce"), type="character", default=NULL,
              help="path to the sce file", metavar="character"),
  make_option(c("-g", "--groupby"), type="character", default=NULL,
              help="name of grouping variable", metavar="character"),
  make_option(c("-a", "--group1"), type="character", default=NULL,
              help="name of cluster 1", metavar="character"),
  make_option(c("-b", "--group2"), type="character", default="All",
              help="name of cluster 2 (default: All)", metavar="character"),
  make_option(c("-k", "--block"), type="character", default="All",
              help="name of blocking factor", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="All",
              help="basic outpath", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opts = parse_args(opt_parser);


# Load data
sce <- readRDS(opts$sce)
sce.edgeR <- prepareEdgeR(sce, opts$groupby, opts$group1,
                          opts$group1, opts$group2, opts$block)


out <- runEdgeR(sce.edgeR,
                min_detection_rate_per_group=0.1,
                max_detection_rate_per_group=1,
                min.logFC=2,
                threshold_fdr=0.1)

saveRDS(out, opts$out)


