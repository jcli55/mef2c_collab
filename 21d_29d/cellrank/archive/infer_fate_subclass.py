import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd

adata = scv.read('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/scvelo/scvelo_result.h5ad')
# Adding subclass annotations to the object
#metadata = pd.read_csv('/storage/singlecell/jeanl/organoid/csv/metadata_clean.csv', index_col = 0)
#list = []
#for barcode in adata.obs_names:
#    list.append(metadata.loc[barcode, 'subclass'])
#adata.obs['subclass'] = list

vk = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.compute_macrostates(n_states=26, cluster_key="subclass")
g.plot_macrostates(which="all", legend_loc="right", s=100, save='macrostates_subclass.png')

g.plot_macrostate_composition(key="subclass", figsize=(7, 4), save='macrostate_composition_class.png')
g.plot_macrostate_composition(key="age", figsize=(7, 4), save='macrostate_composition_age.png')

g.plot_coarse_T(annotate=False, save='coarse_t.png')

g.set_terminal_states(states=['RBC','S_Cone','Rod_1','MG_1','OFF-BC','ML_Cone','Rod_2','MG_2'])
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

g.set_initial_states(states=['PRPC_1','PRPC_2','PRPC_3','PRPC_4','PRPC_5','PRPC_6','PRPC_7','PRPC_8'])
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_subclass.h5ad")