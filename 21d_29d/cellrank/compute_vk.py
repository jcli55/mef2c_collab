import joblib
import cellrank as cr
import scvelo as scv

# CellRank meets RNA velocity

adata = scv.read('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad')
adata.layers['velocity'] = adata.layers['velo_s']

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

joblib.dump(vk, "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/velocity_kernel.h5ad")
