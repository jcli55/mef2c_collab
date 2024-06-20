import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv

# This is part 2 of the multivelo script after running seurat_wnn.R
# Part 1 is in run_multivelo.py
# Adapted for Mef2c collab Jean Li 6/10/2024

# Part 2 ------------------------------------
adata_rna = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna_GC0_normalized_downsampled.h5ad")
adata_atac = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_atac_normalized.h5ad")

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/seurat_wnn/nn_idx.txt", delimiter=',')
nn_dist = np.loadtxt("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/seurat_wnn/nn_dist.txt", delimiter=',')
nn_cells = pd.Index(pd.read_csv("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/seurat_wnn/nn_cells.txt", header=None)[0])

# Make sure cell names match.
#np.all(nn_cells == adata_atac.obs_names)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

# Running Multi-omic Dynamic Model - This will take a while. Parallelization is high recommended.
# mv.settings.VERBOSITY = 0

adata_result = mv.recover_dynamics_chrom(adata_rna, 
                                         adata_atac, 
                                         max_iter=5, 
                                         init_mode="invert",
                                         parallel=True,
                                         save_plot=False,
                                         rna_only=False,
                                         fit=True,
                                         n_anchors=500, 
                                        )

# Plot velocity stream and latent time
scv.pp.neighbors(adata_result)

mv.velocity_graph(adata_result)
mv.latent_time(adata_result)

# Add some metadata labels to the object
age = []
sex = []
geno = []
for sam in adata_result.obs.sampleid:
    age.append('21d') if '21d' in sam else age.append('29d')
    sex.append('f') if '_f_' in sam else sex.append('m')
    geno.append('con') if 'con' in sam else geno.append('ko')
adata_result.obs['age'] = age
adata_result.obs['sex'] = sex
adata_result.obs['geno'] = geno

stage = []
early_cutoff = np.quantile(adata_result.obs.latent_time, 0.05)
late_cutoff = np.quantile(adata_result.obs.latent_time, 0.95)
# The cutoff values are the 0.05 and 0.95 quartile's repectively for latent time
for i in adata_result.obs.latent_time:
    if i < early_cutoff:
        stage.append('early')
    elif i > late_cutoff:
        stage.append('late')
    else:
        stage.append('intermediate')
adata_result.obs['development_stage'] = stage

# Save the result for use later on
adata_result.write("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad")
adata_result.obs.to_csv("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/metadata.csv")

sc.pl.umap(adata_result, color=['sampleid', 'age', 'sex', 'geno', 'development_stage'], frameon=False, ncols=2, size=40, save='_multivelo_metadata.png')
mv.velocity_embedding_stream(adata_result, basis='umap', color='sampleid', save='velocity_stream.png')
scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=80, save='latent_time.png')
