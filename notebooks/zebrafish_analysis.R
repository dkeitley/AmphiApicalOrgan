
# Zebrafish analysis

# Load data
zeb.sce <- loadWagner2018(wagner.path)

AO_genes <- c("foxq2","six3a","six3b","six6a","six6b",
                  "fzd5","fzd8a","fzd8b","sfrp1a","sfrp1b",
                  "sfrp2","sfrp5","dkk3a","dkk3b","dkk1a",
                  "dkk1b","fezf1","fezf2","rx1","rx2","rx3",
                  "nkx2.1","otx2","otx1a","otx1b","otx5","lhx2a","lhx2b",
                  "lhx9")

other.genes <- c("otpa","otpa","wnt8a","wnt8b","wnt11","bsx","gsc")



plotCelltypeHeatmap(zeb.sce,group_by="celltype",genes=apical.genes)

zeb.sce$celltype[zeb.sce$celltype=="NaN"] <- "Other"

zeb.sce$celltype_stage <- paste0(zeb.sce$celltype,"_",zeb.sce$stage)

plotCelltypeStageHeatmap(zeb.sce,group_by = "celltype_stage",genes=apical.genes,
                         celltype_order=c("Forebrain / Optic","Midbrain",
                                          "Hindbrain / Spinal Cord","Neural Crest",
                                          "Epidermal","Mesoderm","Endoderm",
                                          "Germline","Pluripotent","Other"),
                         stage_order=unique(zeb.sce$stage))

zeb.sub <- zeb.sce[,zeb.sce$ClusterName %in% c("24hpf-neural - diencephalon ",
                                               "24hpf-neural - diencephalon posterior",
                                               "18hpf-neural - diencephalon ",
                                               "14hpf-neural - diencephalon ",
                                               "10hpf-neural - diencephalon ",
                                               "10hpf-neural - anterior") ]





# Create forebrain heatmap

# Load STITCH graph nodes meta
stitch_meta<- read.csv("G:/My Drive/Postgrad/PhD/Projects/data/Wagner_2018/STITCH/export_directory/G_nodes.csv")
zeb.sce <- reorderSTITCH(zeb.sce,stitch_meta)


# Load AO annotations
AO_annotation <- read.csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/figures/zebrafish/annotations/AO_annotation.csv" ,
                          col.names = c("sitich_index","AO_annotation"))

# Convert from python to R indexing
AO_annotation$sitich_index = AO_annotation$sitich_index + 1

# Add AO annotations
zeb.sce$AO_annotation <- zeb.sce$celltype
zeb.sce$AO_annotation[zeb.sce$stitch_index!=-1] <- AO_annotation[zeb.sce$stitch_index[zeb.sce$stitch_index!=-1],"AO_annotation"]


# Organise heatmap annotations
zeb.sce$heatmap_annotation <- zeb.sce$celltype
zeb.sce$heatmap_annotation[zeb.sce$celltype == "Forebrain / Optic"] <- "FB/OP"
zeb.sce$heatmap_annotation[zeb.sce$celltype == "Midbrain"] <- "MB"
zeb.sce$heatmap_annotation[zeb.sce$celltype == "Hindbrain / Spinal Cord"] <- "HB/SC"
zeb.sce$heatmap_annotation[zeb.sce$celltype == "Neural Crest"] <- "NC"
zeb.sce$heatmap_annotation[zeb.sce$celltype %in% c("Germline","Other",
                                                   "NaN","Pluripotent")] <- "Other"
zeb.sce$heatmap_annotation[is.na(zeb.sce$celltype)] <- "Other"


# Enforce specific cell type ordering
zeb.sce$heatmap_annotation <- factor(zeb.sce$heatmap_annotation,
                                           levels=c("FB/OP","MB",
                                                    "HB/SC","NC",
                                                    "Epidermal",
                                                    "Mesoderm","Endoderm",
                                                    "Other"))

zeb.sce$stage <- factor(zeb.sce$stage,levels=c("4hpf","6hpf","8hpf","10hpf","14hpf",
                                           "18hpf","24hpf"))

# Plot AO gene expression across tissues
plotDotPlotHeatmap(zeb.sce,zeb.sce$heatmap_annotation,zeb.sce$stage,AO_genes,AO_genes,
                               normalisation="z-score",
                               plot_stage_legend=F)



# Subset forebrain clusters
zeb.fb <-zeb.sce[,zeb.sce$AO_annotation %in% c("Telencephalon","Diencephalon",
                                               "Hypothalamus","Differentiating neurons",
                                               "Optic","Floor plate")]

# Plot AO gene expression across forebrain clusters
plotGroupedHeatmap(zeb.fb,"AO_annotation",AO_genes,AO_genes,
                   c("Telencephalon","Optic","Hypothalamus","Diencephalon",
                     "Differentiating neurons"),keep_gene_order = T)





# Differential expression -------------------------------------------------


# ZF – Hyp vs Ectodermal tissues at 14hpf onwards.
celltypes <- c("Hindbrain / Spinal Cord","Forebrain / Optic","Midbrain",
               "Telencephalon","Diencephalon","Hypothalamus","Floor plate",
               "Differentiating neurons","Optic","Neural crest","Epidermal")
zeb.sub <- zeb.sce[,zeb.sce$AO_annotation%in% celltypes & 
                     zeb.sce$stage %in% c("14hpf","18hpf","24hpf")]
markers <- findMarkers(zeb.sub,groups=zeb.sub$AO_annotation,pval.type="all",
                       block=zeb.sub$library_id)
write.table(markers$Hypothalamus,"G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/gene_lists/differential_expression/zebrafish/hypo_vs_ectoderm/hypo_degs_wblock.tsv",
            sep="\t",quote = F)


# ZF – Second Prosencephalon (Hyp + Retina + Telencephalon) vs Ectodermal Tissues 
# at 14hpf onwards
zeb.sub$DE_annotation <- zeb.sub$AO_annotation
zeb.sub$DE_annotation[zeb.sub$DE_annotation %in% c("Hypothalamus",
                                                   "Optic","Telencephalon")] <- "Secondary Prosencephalon"
markers <- findMarkers(zeb.sub,groups=zeb.sub$DE_annotation,pval.type="all",
                       block=zeb.sub$library_id)
write.table(markers$`Secondary Prosencephalon`,"G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/gene_lists/differential_expression/zebrafish/prosencephalon_vs_ectoderm/degs_wblock.tsv",
            sep="\t",quote = F)


# ZF – Hyp vs Nervous
NS_celltypes <- c("Hindbrain / Spinal Cord","Forebrain / Optic","Midbrain",
               "Telencephalon","Diencephalon","Hypothalamus","Floor plate",
               "Differentiating neurons","Optic","Neural crest")
zeb.sub <- zeb.sce[,zeb.sce$AO_annotation%in% NS_celltypes & 
                     zeb.sce$stage %in% c("14hpf","18hpf","24hpf")]
markers <- findMarkers(zeb.sub,groups=zeb.sub$AO_annotation,pval.type="all",
                       block=zeb.sub$library_id)
write.table(markers$Hypothalamus,"G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/gene_lists/differential_expression/zebrafish/hypo_vs_nervous/degs_wblock.tsv",
            sep="\t",quote = F)


#  ZF – Second Prosencephalon (Hyp + Retina + Telencephalon) vs Nervous

NS_celltypes <- c("Hindbrain / Spinal Cord","Forebrain / Optic","Midbrain",
                  "Telencephalon","Diencephalon","Hypothalamus","Floor plate",
                  "Differentiating neurons","Optic","Neural crest")
zeb.sub <- zeb.sce[,zeb.sce$AO_annotation%in% NS_celltypes & 
                     zeb.sce$stage %in% c("14hpf","18hpf","24hpf")]
zeb.sub$DE_annotation <- zeb.sub$AO_annotation
zeb.sub$DE_annotation[zeb.sub$DE_annotation %in% c("Hypothalamus",
                                                   "Optic","Telencephalon")] <- "Secondary Prosencephalon"
markers <- findMarkers(zeb.sub,groups=zeb.sub$DE_annotation,pval.type="all",
                       block=zeb.sub$library_id)
write.table(markers$`Secondary Prosencephalon`,"G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/gene_lists/differential_expression/zebrafish/prosencephalon_vs_nervous/degs_wblock.tsv",
            sep="\t",quote = F)










# Get hypothalamus markers
markers <- findMarkers(zeb.sce,groups=zeb.sce$AO_annotation,pval.type="some",
                       block=zeb.sce$library_id)
hypo_markers <- markers$Hypothalamus

# Get Forebrain markers
zeb.sce$AO_annotation_fb <- zeb.sce$AO_annotation
zeb.sce$AO_annotation_fb[zeb.sce$AO_annotation_fb %in% c("Telencephalon","Hypothalamus","Optic",
                                                         "Differentiating neurons","Floor plate")] <- "Forebrain"
fb_markers <- findMarkers(zeb.sce,groups=zeb.sce$AO_annotation_fb,pval.type="some",
                       block=zeb.sce$library_id)
fb_markers <- fb_markers$Forebrain

write.table(hypo_markers,"G:/My Drive/Postgrad/PhD/Projects/apical_organ/gene_lists/differential_expression/zebrafish/some/hypo_degs_wblock.tsv",
            sep="\t",quote = F)

zeb.sce$AO_annotation[zeb.sce$AO_annotation %in% c("Hypothalamus","Optic")] <- "Hypo/Optic"
hypo_optic_markers <- findMarkers(zeb.sce,groups=zeb.sce$AO_annotation,pval.type="all")



# Plot volcano plot
library(EnhancedVolcano)
EnhancedVolcano(hypo_markers,x="summary.logFC",y="FDR",
                ylab="-log2 FDR",
                xlab="log2 FC ",
                lab=rownames(hypo_markers),
                selectLab = rownames(hypo_markers)[1:10],
                labSize = 4,
                drawConnectors = T,arrowheads = F,
                title="",subtitle = "",caption = "",
                legendPosition = "none",
                shape=16) 
