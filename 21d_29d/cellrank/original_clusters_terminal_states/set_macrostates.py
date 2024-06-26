import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd
import numpy as np

adata = scv.read('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad')

vk = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/velocity_kernel.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

adata.obs['initialstates'] = adata.obs.development_stage.replace(
	{
		"late": np.nan,
		"intermediate": np.nan,
	}
)

g.set_initial_states(states = adata.obs.initialstates)
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial_orig.png')

adata.obs['terminalstates'] = adata.obs.orig_stages.replace(
	{
		'intermediate1': np.nan, 
        'intermediate0': np.nan, 
        'early0': np.nan, 
        'early2': np.nan, 
        'early1': np.nan, 
        'late5': np.nan, 
        'intermediate5': np.nan, 
        'early5': np.nan
	}
)

g.set_terminal_states(states = adata.obs.terminalstates)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal_orig.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont_orig.png')

joblib.dump(g, "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator_orig_stage.h5ad")
