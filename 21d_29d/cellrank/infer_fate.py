import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd

adata = scv.read('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad')

vk = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/velocity_kernel.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.fit(n_states=[4, 12])
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100, save='macrostates.png')

g.plot_macrostate_composition(key='age', save='macrostate_age_composition.png')
g.plot_macrostate_composition(key='sex', save='macrostate_sex_composition.png')
g.plot_macrostate_composition(key='geno', save='macrostate_geno_composition.png')

g.predict_terminal_states(allow_overlap=True)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')
g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

g.predict_initial_states(allow_overlap=True)
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

adata.obs["macrostates"] = g.macrostates
adata.uns["macrostates_colors"] = g.macrostates_memberships.colors
adata.write('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/multivelo_result.h5ad')
adata.obs.to_csv('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/metadata.csv')

joblib.dump(g, "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator.h5ad")
