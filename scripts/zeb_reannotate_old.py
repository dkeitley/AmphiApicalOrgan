# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 07:50:19 2021

@author: Daniel Keitley
"""

import scanpy as sc
from scipy.sparse import dok_matrix
import scipy.sparse as ss
import pandas as pd
import numpy as np

sc.set_figure_params(dpi_save=300)
sc.settings.figdir = "G:/My Drive/Postgrad/PhD/Projects/ApicalOrgan2021/plots/"


stitch_path = "G:/My Drive/Postgrad/PhD/Projects/data/Wagner_2018/STITCH/export_directory/"

zeb = sc.read_csv(stitch_path + "counts.csv")
zeb = sc.read_h5ad(stitch_path + "zeb_temp.h5ad")

# Add gene names
genes = pd.read_csv(stitch_path + 'genes.txt', delimiter = "\t",
                        header=None,dtype='category')
zeb.var.index = genes[0]

# Add cell type labels
celltypes = pd.read_csv(stitch_path + 'cell_IDs_names.txt', delimiter = "\t",
                        header=None,dtype='category')
zeb.obs["celltype"] = celltypes[0].values

tissues = pd.read_csv(stitch_path + 'branch_IDs_names.txt', delimiter = "\t",
                        header=None,dtype='category')
zeb.obs["tissue"] = tissues[0].values


timepoints = pd.read_csv(stitch_path + 'timepoints.txt',
                         header=None,dtype='category')
timepoints = timepoints[0].cat.rename_categories(["4hpf","6hpf","8hpf","10hpf","14hpf","18hpf","24hpf"])

zeb.obs["stage"] = timepoints.values

# Add cell type colours
#colours = pd/

# Add graph from STITCH
edges = pd.read_csv(stitch_path + "edges.csv",sep=";",header=None)

adj_mat = ss.dok_matrix((zeb.n_obs,zeb.n_obs))
for index, row in edges.iterrows():
    adj_mat[row[0],row[1]] = 1

zeb.uns['neighbors'] = {}
neighbors_dict = zeb.uns['neighbors']
neighbors_dict['connectivities_key'] = 'connectivities'
zeb.obsp["connectivities"] = adj_mat.tocsr()
zeb.uns["neighbors"] = neighbors_dict

sc.tl.draw_graph(zeb)
sc.pl.draw_graph(zeb,color="tissue")


zeb_fb = zeb[zeb.obs["tissue"].isin(["Forebrain / Optic"]),]
             
zeb_dc = zeb[zeb.obs["celltype"].isin(["24hpf-neural - diencephalon ",
                                         "24hpf-neural - diencephalon posterior",
                                         "18hpf-neural - diencephalon ",
                                         "14hpf-neural - diencephalon "]),]

zeb_dc = zeb_fb[zeb_fb.obs["celltype"].isin(["24hpf-neural - diencephalon ",
                                         "18hpf-neural - diencephalon " ]),]

celltype_colours = {"24hpf-neural - diencephalon ":"blue",
  "24hpf-neural - diencephalon posterior ":"red",
  "18hpf-neural - diencephalon ":"green",
  "14hpf-neural - diencephalon ":"yellow",
  "10hpf-neural - anterior":"purple",
  "10hpf-neural - diencephalon ":"black"
 }

zeb_dc.uns["celltype_colors"][2] = '#D5CE81'
zeb_dc.uns["celltype_colors"][0] = '#A4D8CA'
zeb_dc.uns["celltype_colors"][1] = '#93502D'
zeb_dc.uns["celltype_colors"][3] = 'green'
zeb_dc.uns["celltype_colors"][4] = '#5394BF'

sc.tl.draw_graph(zeb_dc)
sc.pl.draw_graph(zeb_dc,color='celltype',size=10)


zeb_dc = zeb[zeb.obs["celltype"].isin(["24hpf-neural - diencephalon ",
                                         "18hpf-neural - diencephalon ",
                                         "14hpf-neural - diencephalon "]),]


zeb_dc = sc.read_h5ad("G:/My Drive/Postgrad/PhD/Projects/apical_organ/zeb reannotation - first attempt/zeb_dc.h5ad")


sc.tl.leiden(zeb_dc)

sc.pl.draw_graph(zeb_dc,color=['celltype','leiden','barhl2','fezf2','nkx2.1',
                               'pitx2','six3a','rx3','foxb1a','dbx1a','pitx3',
                               "irx3a"],size=50,
                 legend_loc='on data')


zeb_dc.obs["leiden"][~zeb_dc.obs["leiden"].isin(["0","1","2","3","4","5","6","7","8"])] = "9"

zeb_dc.obs["leiden"][~zeb_dc.obs["leiden"].isin(["0","1","2","3","4","5"])] = "6"
zeb_dc.obs["leiden"] = zeb_dc.obs["leiden"].cat.remove_unused_categories()
dc_filt = zeb_dc[~zeb_dc.obs["leiden"].isin(["9"])]


sc.tl.draw_graph(dc_filt)

sc.pl.draw_graph(dc_filt,color=['celltype','leiden','barhl2','fezf2','nkx2.1',
                               'pitx2','six3a','rx3','foxb1a','dbx1a','pitx3',
                               "irx3a"],size=50,
                 legend_loc='on data')

dc_filt.obs["edited_celltype"] = dc_filt.obs["celltype"].astype(str)
dc_filt.obs["edited_celltype"][dc_filt.obs["leiden"].isin(["0","1","2","4","5"])] = "Hypothalamus"
sc.pl.draw_graph(dc_filt,color=['celltype','edited_celltype'],size=50)

zeb.obs["edited_celltype"] = 
zeb.obs["edited_celltype"][dc_filt.obs["leiden"].isin(["0","1","2","4","5"]).index] = "Hypothalamus"

# Save hypothalamus annotation
np.savetxt("G:/My Drive/Postgrad/PhD/Projects/apical_organ/figures/zebrafish/annotations/hypo_annotation_stitch_index.csv",
           dc_filt.obs[dc_filt.obs["edited_celltype"]=="Hypothalamus"].index.values.astype(int),
           fmt='%i')


# Make temporary dot plot
late_fb = zeb_fb[zeb_fb.obs["stage"].isin(["24hpf","18hpf","14hpf"]),]

# Reannotate
late_fb.obs["edited_celltype"] = ""
late_fb.obs["edited_celltype"][late_fb.obs["celltype"].isin(['14hpf-optic primordium',
                                                             '18hpf-optic cup',
                                                             '18hpf-retina pigmented epithelium',
                                                             '24hpf-optic cup',
                                                             '24hpf-retina pigmented epithelium'])] = "Optic"

late_fb.obs["edited_celltype"][late_fb.obs["celltype"].isin(['14hpf-differentiating neurons - rohon beard',
                                                             '18hpf-differentiating neurons - dlx',
       '18hpf-differentiating neurons - eomesa',
       '18hpf-differentiating neurons - rohon beard',
       '24hpf-differentiating neurons - dlx',
       '24hpf-differentiating neurons - eomesa',
       '24hpf-differentiating neurons - hmx',
       '24hpf-differentiating neurons - rohon beard'
       ])] = "Differentiating neurons"

late_fb.obs["edited_celltype"][late_fb.obs["celltype"].isin(['14hpf-neural -  floorplate',
                                                             '18hpf-neural -  floorplate',
                                                             '24hpf-neural - floorplate'])] = "Floor plate"

late_fb.obs["edited_celltype"][late_fb.obs["celltype"].isin(['14hpf-neural - telencephalon',
                                                             '18hpf-neural - telencephalon',
                                                             '24hpf-neural - telencephalon'])] = "Telencephalon"
                               
late_fb.obs["edited_celltype"][late_fb.obs["celltype"].isin(['14hpf-neural - diencephalon ',
                                                             '18hpf-epiphysis',
                                                             '18hpf-neural - diencephalon ',
                                                             '24hpf-neural - diencephalon ',
       '24hpf-neural - diencephalon posterior'])] = "Diencephalon"



late_fb.obs["edited_celltype"][dc_filt.obs[dc_filt.obs["edited_celltype"]=="Hypothalamus"].index] = "Hypothalamus"
                

apical_genes = ["foxq2","six3a","six3b","six6a","six6b",
                  "fzd5","fzd8a","fzd8b","sfrp1a","sfrp1b",
                  "sfrp2","sfrp5","dkk3a","dkk3b","dkk1a",
                  "dkk1b","fezf1","fezf2","rx1","rx2","rx3",
                  "nkx2.1","otx2","otx1a","otx1b","otx5","otpa",
                  "otpa","lhx2a","lhx2b",
                  "lhx9"]
sc.pl.dotplot(late_fb,apical_genes,groupby="edited_celltype",
              swap_axes=False,
              categories_order=["Telencephalon","Optic","Hypothalamus",
                                "Diencephalon","Differentiating neurons","Floor plate"],
              save="_zeb_forebrain_dotplot_temp2.pdf")


sc.pl.draw_graph(late_fb,color="edited_celltype",
                 save="_zeb_late_fb_fa.pdf")







sc.pl.violin(dc_filt,['barhl2','fezf2','nkx2.1',
                     'pitx2','six3a','rx3','foxb1a',
                     'dbx1a'],"leiden")


sc.tl.rank_genes_groups(dc_filt, 'leiden')
sc.pl.rank_genes_groups(dc_filt, n_genes=25, sharey=False)


zeb_fb.obs["dc_annotation"] = ""
zeb_fb.uns["celltype_colors"]["24hpf-neural - diencephalon posterior"] = '#5394BF'

possible_hyp_cells = zeb_dc.obs["leiden"]
zeb_fb.obs["dc_annotation"][zeb_dc.obs.index] = zeb_dc.obs["leiden"]

zeb_dc.write_h5ad("G:\\My Drive\\Postgrad\\PhD\\Projects\\apical_organ\\zeb reannotation - first attempt\\zeb_dc.h5ad")




# apical_genes = ["foxq2","six3a","six3b","six6a","six6b",
#                   "fzd5","fzd8a","fzd8b","sfrp1a","sfrp1b",
#                   "sfrp2","sfrp5","dkk3a","dkk3b","dkk1a",
#                   "dkk1b","fezf1","fezf2","rx1","rx2","rx3",
#                   "nkx2.1","otx2","otx1a","otx1b","otx5","otpa",
#                   "otpa","lhx2a","lhx2b"]

# sc.pl.dotplot(zeb,apical_genes,groupby="tissue",swap_axes=True,
#               cmap="viridis")

import seaborn as sns
