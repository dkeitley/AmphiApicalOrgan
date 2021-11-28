# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 12:08:08 2021

@author: Daniel Keitley
"""

import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pathlib import Path

import plotly.io as pio
pio.renderers.default='browser'

sc.set_figure_params(dpi=300,dpi_save=300)

import os
os.chdir("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/")

from scripts.annotation_utils import *


# Load data
echino = sc.read_h5ad("G:/My Drive/Postgrad/PhD/Projects/data/Foster_2020/echino_orig.h5ad")
new_pca = pd.read_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/data-out/batch_correction/echino_fastMNN_pcs.tsv",sep="\t")

echino.obsm["X_pca"] = np.array(new_pca)

AO_genes = pd.read_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/data-in/echino_heatmap_genes.tsv",sep=";")
AO_genes = AO_genes.append({"Gene":"TPH","Code":"LOC578903"},ignore_index=True)

# Add gene names to AO genes
echino.var["gene_name"] = ""
echino.var["gene_name"] = echino.var["gene_name"] .astype(str)
echino.var["gene_name"].loc[echino.var.index.isin(AO_genes["Code"])] = AO_genes["Gene"].values


# Compute UMAP
sc.pp.neighbors(echino, n_neighbors=75, use_rep="X_pca", n_pcs=50)
sc.tl.umap(echino, min_dist=0.9)


# Export UMAP
umap_df = pd.DataFrame(echino.obsm["X_umap"], columns=["UMAP_1","UMAP_2"])
umap_df.to_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/umap.tsv", sep="\t", index=False)


# Load UMAP
umap_df = pd.read_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/umap.tsv", sep="\t")
echino.obsm["X_umap"] = np.array(umap_df)


# Visualise existing annotations and expression
sc.pl.umap(echino, color="AO_annotation")
sc.pl.umap(echino, color="AO_annotation", groups=["Serotonergic neurons", "Late AP",
                                                  "Early AP", "Other neural"])

sc.pl.umap(echino, color=AO_genes["Code"],title=AO_genes["Gene"],color_map="viridis_r",size=2)


# Perform clustering
sc.tl.leiden(echino, resolution=1, key_added="leiden_res1",random_state=0)
sc.tl.leiden(echino, resolution=2, key_added="leiden_res2",random_state=0)
sc.tl.leiden(echino, resolution=3, key_added="leiden_res3",random_state=0)
sc.tl.leiden(echino, resolution=5, key_added="leiden_res5",random_state=0)
sc.tl.leiden(echino, resolution=8, key_added="leiden_res8",random_state=0)


# Make annotation plots
makeAnnotationPlots(echino, 
                    clusters = ["leiden_res1", "leiden_res2","leiden_res3","leiden_res5","leiden_res8"],
                    model_predictions = "AO_annotation", other_obs=["orig.ident"],
                   export_dir = "G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/")



# Subset apical organ region
echino_ap =  echino[echino.obs["leiden_res2"].isin(['17','7','9','10','11','29','18','24','2']),]
sc.tl.leiden(echino_ap,resolution=2,key_added="ap_leiden_res2",random_state=0)
sc.tl.leiden(echino_ap,resolution=3,key_added="ap_leiden_res3",random_state=0)
sc.tl.leiden(echino_ap,resolution=4,key_added="ap_leiden_res4",random_state=0)
sc.tl.leiden(echino_ap,resolution=5,key_added="ap_leiden_res5",random_state=0)


# Compute 3D UMAP
umap_3d = sc.pp.neighbors(echino_ap,n_neighbors=25, use_rep="X_pca", random_state=0, copy=True)
umap_3d = sc.tl.umap(umap_3d, min_dist=0.5, n_components=3, random_state=0, copy=True)
umap_3d.obs["ap_leiden_res2"] = echino_ap.obs["ap_leiden_res2"]
umap_3d.obs["ap_leiden_res3"] = echino_ap.obs["ap_leiden_res3"]
umap_3d.obs["ap_leiden_res4"] = echino_ap.obs["ap_leiden_res4"]
umap_3d.obs["ap_leiden_res5"] = echino_ap.obs["ap_leiden_res5"]


makeAnnotationPlots(echino_ap, 
                    clusters = ["ap_leiden_res4"],
                    model_predictions = "AO_annotation", other_obs=["orig.ident"],
                   export_dir = "G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/echino_ap/")


# Visualise in 3D
plotObsUMAP3d(umap_3d,"ap_leiden_res5",
                     hover_obs=["AO_annotation","ap_leiden_res5","orig.ident","leiden_res8"],
                    figsize=(700,600))


# Plot expression across clusters
sc.pl.matrixplot(echino_ap, var_names=AO_genes["Gene"], groupby="ap_leiden_res5",
                 standard_scale='var',gene_symbols="gene_name")


# Other neural
sc.pl.umap(echino, color=["ebr1","egf1","LOC575288","LOC577619","LOC587837",
                          "LOC576524","LOC373478","delta",
                          "LOC577601", "LOC105445683",""
                          "SoxB1","SoxB2"], color_map="viridis_r")


# Assign cell types to clusters
echino.obs["assigned_celltype"] = ""
echino_ap.obs["assigned_celltype"] = ""

# Perform assignment
assignCelltypeAdata(echino, "leiden_res5", ["0","7","36","67"], "assigned_celltype", "Other neural")
assignCelltypeAdataSubset(echino, echino_ap, "ap_leiden_res5", ["37","24","19","17"], "assigned_celltype", "Early AP")
assignCelltypeAdataSubset(echino, echino_ap, "ap_leiden_res5", ["40","9","32","29","8","4"], "assigned_celltype", "Late AP")
assignCelltypeAdataSubset(echino, echino_ap, "ap_leiden_res5", ["39"], "assigned_celltype", "Serotonergic neurons")

# Visualise new assignements
sc.pl.umap(echino_ap, color="assigned_celltype")

umap_3d.obs["assigned_celltype"] = echino_ap.obs["assigned_celltype"]
plotObsUMAP3d(umap_3d,"assigned_celltype",
                     hover_obs=["AO_annotation","ap_leiden_res5","orig.ident"],
                    figsize=(700,600))

sc.pl.umap(echino, color="assigned_celltype", groups=["Serotonergic neurons", "Late AP",
                                                  "Early AP"])


# Export new cell types
echino.obs["assigned_celltype"].to_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/AO_annotations.tsv",sep="\t")

echino.obs["assigned_celltype"] = echino.obs["assigned_celltype"].astype(str)
echino.obs["AO_annotation"].loc[echino.obs["assigned_celltype"]!=""] = echino.obs["assigned_celltype"].loc[echino.obs["assigned_celltype"]!=""]
echino.obs["AO_annotation"].to_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/AO_annotations_v2.tsv",sep="\t")



sc.tl.rank_genes_groups(echino, groupby="leiden_res2", groups=["26"])
sc.pl.rank_genes_groups(echino,groups=["26"])



AO_annotations = pd.read_csv("G:/My Drive/Postgrad/PhD/Projects/data/Foster_2020/Gattoni2021_annotation.tsv",sep="\t",header=None)

echino.obs["AO_annotation"] = AO_annotations.iloc[:,0].values
echino.obs["AO_annotation"][echino.obs["AO_annotation"].isin(["Serotonergic neurons", "Late AP",
                                                  "Early AP", "Other neural"])] = "Other"


echino.obs["AO_annotation"].to_csv("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/echino_annotation/AO_annotations.tsv",sep="\t")
