import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd
import numpy as np

adata = scv.read('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad')

adata.obs["maturityclass"] = adata.obs.subclass.replace(
    {
        "ML_Cone": "Cone",
        "S_Cone": "Cone",
        "OFF-BC": "BC",
        "ON-BC": "BC",
        "RBC": "BC",
        "GABAergic": "AC",
        "Glycinergic": "AC",
        "SACs": "AC",
        "dual ACs": "AC",
        "HC0": "HC",
        "HC1": "HC"
    }
)

vk = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel_multivelo.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.compute_macrostates(n_states=14, cluster_key="maturityclass")
g.plot_macrostates(which="all", legend_loc="right", s=100, save='macrostates.png')

g.plot_macrostate_composition(key="maturityclass", figsize=(7, 4), save='macrostate_composition_maturityclass.png')
g.plot_macrostate_composition(key="age", figsize=(7, 4), save='macrostate_composition_age.png')

g.plot_coarse_T(annotate=False, save='coarse_t.png')

g.set_initial_states(states=['PRPC_1','PRPC_5'])
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

adata.obs['terminalstates'] = adata.obs.maturityclass.replace(
    {
        "Cone Precursor": np.nan,
        "Rod Precursor": np.nan,
        "BC Precursor": np.nan,
        "AC Precursor": np.nan,
        "HC Precursor": np.nan,
        "PRPC": np.nan,
        "NRPC": np.nan
    }
)

g.set_terminal_states(states = adata.obs.terminalstates)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

# write macrostates to AnnData
adata.obs["macrostates"] = g.macrostates
adata.uns["macrostates_colors"] = g.macrostates_memberships.colors

adata.write('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad')

joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo_scvi_states.h5ad")