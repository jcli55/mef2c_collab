import joblib
import cellrank as cr
import scvelo as scv

# CellRank meets RNA velocity

adata = scv.read('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/scvelo/scvelo_result.h5ad')

vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()

joblib.dump(vk, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel.h5ad")
