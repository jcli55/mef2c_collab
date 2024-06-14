import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv

# Runs part 1 of the multivelo pipeline
# Adapted to run on Suyang's data - 6/7/24
# Had to add a few lines to format the obs_names the same for the intersects

# Demo params --------------
# scv.settings.verbosity = 3
# scv.settings.presenter_view = True
# scv.set_figure_params('scvelo')
# pd.set_option('display.max_columns', 100)
# pd.set_option('display.max_rows', 200)
# np.set_printoptions(suppress=True)

# Part 1 -------------------------------------
# Read in and prepare the merged rna and atac data
# RNA
adata_rna = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/merged_spliced_unspliced.h5ad")
adata_rna.var_names_make_unique()
names = []
for name in adata_rna.obs_names: #Make the obs_names match the other objects
    if '21d_f_con' in name: #Check each sampleid in the obs_names and add on the correct one to the barcode
        names.append(f'21d_f_con#{name[len(name) - 20: ]}')
    elif '21d_f_ko' in name:
        names.append(f'21d_f_ko#{name[len(name) - 20: ]}')
    elif '21d_m_con' in name:
        names.append(f'21d_m_con#{name[len(name) - 20: ]}')
    elif '21d_m_ko' in name:
        names.append(f'21d_m_ko#{name[len(name) - 20: ]}')
    elif '29d_m_con' in name:
        names.append(f'29d_m_con#{name[len(name) - 20: ]}')
    elif '29d_m_ko' in name:
        names.append(f'29d_m_ko#{name[len(name) - 20: ]}')
adata_rna.obs_names = names
#names = []
#for name in adata_rna.obs_names:
#    names.append(name[len(name) - 20: ]) # Only take the last 20 characters of each obs_name (the 16bp barcode and the -1-0 tag)
#adata_rna.obs['cellid'] = adata_rna.obs_names
#adata_rna.obs_names = names
#adata_rna.obs_names_make_unique()

# Transfer the cell labels to the spliced + unspliced data
filtered_annotated_rna = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna_GC0.h5ad")
names = []
for name in filtered_annotated_rna.obs_names:
    names.append(f'{name}-0')
filtered_annotated_rna.obs_names = names
# names = []
# for name in filtered_annotated_rna.obs_names:
#     name = f'{name}-0'
#     names.append(name[len(name) - 20: ])
#filtered_annotated_rna.obs_names = names
#filtered_annotated_rna.obs_names_make_unique()

filtered_cells = pd.Index(np.intersect1d(adata_rna.obs_names, filtered_annotated_rna.obs_names))
adata_rna = adata_rna[filtered_cells]

sampleid = []
for name in adata_rna.obs_names:
    sampleid.append(filtered_annotated_rna.obs.loc[name]['sampleid'])
adata_rna.obs['sampleid'] = sampleid
adata_rna.obs['sampleid'] = adata_rna.obs['sampleid'].astype('category')

# ATAC
adata_atac = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_atac.h5ad")
names = []
for name in adata_atac.obs_names:
    names.append(f'{name}-0')
adata_atac.obs_names = names
# for name in adata_atac.obs_names:
#     name = f'{name}-0'
#     names.append(name[len(name) - 20: ])
# adata_atac.obs['cellid'] = adata_atac.obs_names
# adata_atac.obs_names = names
# adata_atac.obs_names_make_unique()

# Get shared cells and genes between rna and atac and subset on just these
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

# Normalization for rna and atac
scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50) # Might need to move to after subsetting bc it will corrupt the neighbors graph

mv.tfidf_norm(adata_atac)

# Subset cells and genes (to reduce runtime and memory requirements)
sc.pp.highly_variable_genes(adata_rna, flavor="seurat_v3", n_top_genes=2000, subset=True)
sc.pp.subsample(adata_rna, n_obs=10000)

# Another round of subset to get atac to match
shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

# Transfer the umap coordinates to the spliced/unspliced data
# Create a dictionary of sampleid/barcodes mapped to their umap coordinate
dict = {}
for i in range(0, len(filtered_annotated_rna)):
    dict[filtered_annotated_rna.obs_names[i]] = filtered_annotated_rna.obsm['X_umap'][i]

# Recreate the X_umap stacked array using the order the velocyto data was in
array = dict[adata_rna.obs_names[0]]
for i in range(1, len(adata_rna)):
    array = np.vstack([array, dict[adata_rna.obs_names[i]]])

# Set the umap
adata_rna.obsm['X_umap'] = array

adata_rna.write("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna_GC0_normalized_downsampled.h5ad")
adata_atac.write("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_atac_normalized.h5ad")

# Generate the filtered cells txt for the R script
adata_rna.obs_names.to_frame().to_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/seurat_wnn/filtered_cells.txt', header=False, index=False)

# -------------------------------------------
# Run Seurat WNN script to generate WNN files
# -------------------------------------------

# # Moved part 2 to its own script run_multivelo_pt2.py
# # Part 2 ------------------------------------
# adata_rna = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_rna_normalized_downsampled.h5ad")
# adata_atac = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/adata_atac_normalized_downsampled.h5ad")

# # Read in Seurat WNN neighbors.
# nn_idx = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_idx.txt", delimiter=',')
# nn_dist = np.loadtxt("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_dist.txt", delimiter=',')
# nn_cells = pd.Index(pd.read_csv("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/seurat_wnn/nn_cells.txt", header=None)[0])

# # Make sure cell names match.
# #np.all(nn_cells == adata_atac.obs_names)

# mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

# # Running Multi-omic Dynamic Model - This will take a while. Parallelization is high recommended.
# # mv.settings.VERBOSITY = 0

# adata_result = mv.recover_dynamics_chrom(adata_rna, 
#                                          adata_atac, 
#                                          max_iter=5, 
#                                          init_mode="invert",
#                                          parallel=True,
#                                          save_plot=False,
#                                          rna_only=False,
#                                          fit=True,
#                                          n_anchors=500, 
#                                         )

# # Save the result for use later on
# adata_result.write("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad")
