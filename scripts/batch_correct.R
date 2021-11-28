

library(scran)
library(batchelor)

source("../utils/load.R")

set.seed(42)

# Load sea urchin data
echino.sce <- loadFoster2020(foster.path)


# Compute highly variable genes
gene_var <- modelGeneVar(echino.sce)
hvgs <- getTopHVGs(gene_var,n=2000)


out <- fastMNN(echino.sce,batch=echino.sce$stage,k=20,subset.row = hvgs,
               merge.order = rev(levels(echino.sce$stage)))


saveRDS(out,"data-out/batch_correction/echino_fastMNN_out.rds")
write.table(reducedDim(out,"corrected"),"data-out/batch_correction/echino_fastMNN_pcs.tsv",sep="\t",quote=F)
