import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd

adata = scv.read('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/scvelo/scvelo_result.h5ad')
# Adding age annotations to the object
#metadata = pd.read_csv('/storage/singlecell/jeanl/organoid/csv/metadata_clean.csv', index_col = 0)
#list = []
#for barcode in adata.obs_names:
#    list.append(metadata.loc[barcode, 'age'])
#adata.obs['age'] = list

vk = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.compute_macrostates(n_states=12, cluster_key="class")
g.plot_macrostates(which="all", legend_loc="right", s=100, save='macrostates.png')

g.plot_macrostate_composition(key="class", figsize=(7, 4), save='macrostate_composition_class.png')
g.plot_macrostate_composition(key="age", figsize=(7, 4), save='macrostate_composition_age.png')

g.plot_coarse_T(annotate=False, save='coarse_t.png')

g.set_terminal_states(states=['MG','Cone_1','Cone_2','Rod_1','Rod_2','AC','BC'])
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

g.set_initial_states(states=['PRPC_1','PRPC_2'])
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_advanced.h5ad")