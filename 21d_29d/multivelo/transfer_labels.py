import scanpy as sc
import pandas as pd
import numpy as np

# Load the data and metadata to transfer
adata = sc.read_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad")

reintegrated = pd.read_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/GC0.re-integrate.re-subcluster.meta.csv')
metadata = pd.read_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/GC0_metadata.csv')

# Match the formatting of the names of the metadata to the adata and set it as the index
index = []
for name in reintegrated['Unnamed: 0']:
    index.append(name + '-0')

reintegrated['barcode'] = index
metadata['barcode'] = index

reintegrated = reintegrated.set_index('barcode')
metadata = metadata.set_index('barcode')

# Transfer the clusters from the metadata to the adata
reintegrated_clusters = []
orig_clusters = []
for name in adata.obs_names:
    reintegrated_clusters.append(reintegrated.loc[name]['RNA_snn_res.0.1'])
    orig_clusters.append(metadata.loc[name]['RNA_snn_res.0.1'])

adata.obs['reintegrated_clusters'] = reintegrated_clusters
adata.obs['orig_clusters'] = orig_clusters

# Create new stage labels to be used in Cellrank (combines the development stage based on latent time to the clusters that were just transferred)
reintegrated_stages = []
orig_stages = []
for name in adata.obs_names:
    reintegrated_stages.append(adata.obs.loc[name]['development_stage'] + str(adata.obs.loc[name]['reintegrated_clusters']))
    orig_stages.append(adata.obs.loc[name]['development_stage'] + str(adata.obs.loc[name]['orig_clusters']))

adata.obs['reintegrated_stages'] = reintegrated_stages
adata.obs['orig_stages'] = orig_stages

# Create columns for the Umap coordinates and Seurat formatted cell id's to transfer back to the Seurat object
x = []
y = []
for i in range(0, len(adata)):
    x.append(adata.obsm['X_umap'][i][0])
    y.append(adata.obsm['X_umap'][i][1])
adata.obs['umap_x'] = x
adata.obs['umap_y'] = y

list = []
for name in adata.obs_names:
    list.append(i.replace('-1-0', '-1'))
adata.obs['Seurat_id'] = list

adata.obs.to_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/metadata.csv')
adata.write('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad')