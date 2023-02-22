# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:25:13 2021

@author: Daniel Keitley
"""


# Setup
import scanpy as sc
from scipy.sparse import dok_matrix
import scipy.sparse as ss
import pandas as pd
import numpy as np


stitch_path = "G:/My Drive/Postgrad/PhD/Projects/data/Wagner_2018/STITCH/export_directory/"


# Load stitch count data
zeb_stitch = sc.read_csv(stitch_path + "counts.csv")

# Load stich graph
edges = pd.read_csv(stitch_path + "edges.csv",sep=";",header=None)


# Add metadata

# Add genes
genes = pd.read_csv(stitch_path + 'genes.txt', delimiter = "\t",
                        header=None,dtype='category',names=["gene_name"])
zeb_stitch.var = genes["gene_name"]

# Add cell type labels
celltypes = pd.read_csv(stitch_path + 'cell_IDs_names.txt', delimiter = "\t",
                        header=None,dtype='category')
zeb_stitch.obs["celltype"] = celltypes[0].values

tissues = pd.read_csv(stitch_path + 'branch_IDs_names.txt', delimiter = "\t",
                        header=None,dtype='category')
zeb_stitch.obs["tissue"] = tissues[0].values

# Add stage
timepoints = pd.read_csv(stitch_path + 'timepoints.txt',
                         header=None,dtype='category')
timepoints = timepoints[0].cat.rename_categories(["4hpf","6hpf","8hpf","10hpf","14hpf","18hpf","24hpf"])
zeb_stitch.obs["stage"] = timepoints.values




# Create adjacency matrix from edges data frame
adj_mat = ss.dok_matrix((zeb_stitch.n_obs,zeb_stitch.n_obs))
for index, row in edges.iterrows():
    adj_mat[row[0],row[1]] = 1
    

# Add STITCH graph to anndata object
zeb_stitch.uns['neighbors'] = {}
neighbors_dict = zeb_stitch.uns['neighbors']
neighbors_dict['connectivities_key'] = 'connectivities'
zeb_stitch.obsp["connectivities"] = adj_mat.tocsr()
zeb_stitch.uns["neighbors"] = neighbors_dict


# Export as anndata object
zeb_stitch.write("G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/data-out/zeb_stitch.h5ad")


