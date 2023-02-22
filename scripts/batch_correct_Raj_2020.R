suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(batchelor))

comb.sce <- readRDS("comb_sce.rds")

comb.sce <- multiBatchNorm(comb.sce, batch=comb.sce$orig.ident,
                           preserve.single=TRUE)

gene_var <- modelGeneVar(comb.sce, block=comb.sce$orig.ident)
hvgs <- getTopHVGs(gene_var,n=2000)

df = colData(comb.sce)[, c("stage", "orig.ident")]
df$ncells = sapply(df$orig.ident, function(x) sum(df$orig.ident == x))

df$stage = factor(df$stage,levels = unique(zeb.clusts$Stage.time))
df$orig.ident <- as.character(df$orig.ident)
df = df[order(df$stage, df$ncells, decreasing = TRUE),]

merge_order <- rev(lapply(split(df,df$stage), function(x) list(unique(x$orig.ident))))
names(merge_order)<- NULL

out <- fastMNN(comb.sce,batch=comb.sce$orig.ident,k=20,subset.row = hvgs,
               merge.order = merge_order)

saveRDS(out,"batch_correction.rds")
write.table(reducedDim(out,"corrected"),"corrected_pcs.tsv",sep="\t",quote=F)
