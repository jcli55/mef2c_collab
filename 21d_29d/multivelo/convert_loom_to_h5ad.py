import scvelo as scv
import anndata

ldata_file = "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/merged.loom"
output_file = "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/merged_spliced_unspliced.h5ad"

ldata = anndata.read_loom(ldata_file)

ldata.obs.index = [
    x.split(":", 1)[0] + "_" + x.split(":", 1)[1][:-1] + "-1-0" for x in ldata.obs.index
]
index = list(ldata.obs.index)

ldata.obs.index = index
ldata.raw = ldata

ldata.write(output_file)
