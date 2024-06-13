import scanpy as sc
import numpy as np
import pandas as pd

print("Finished importing, starting reading...")
adata = sc.read_h5ad('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna.h5ad')
metadata = pd.read_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/GC0_metadata.csv')
metadata = metadata.set_index('Unnamed: 0')
print("Finished reading files")

intersect = pd.Index(np.intersect1d(adata.obs_names, metadata.index))
adata = adata[intersect]

sampleid = []
age = []
sex = []
geno = []

print("Starting transfer...")
for i in range(0, len(adata)):
    barcode = adata.obs_names[i]
    sampleid.append(metadata.loc[barcode]['orig.ident'])
    age.append(metadata.loc[barcode]['age'])
    sex.append(metadata.loc[barcode]['sex'])
    geno.append(metadata.loc[barcode]['geno'])

adata.obs['sampleid'] = sampleid
adata.obs['age'] = age
adata.obs['sex'] = sex
adata.obs['geno'] = geno
print("Finished transfer")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sampleid")
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

print("Writing file...")
adata.write_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna_GC0.h5ad")
