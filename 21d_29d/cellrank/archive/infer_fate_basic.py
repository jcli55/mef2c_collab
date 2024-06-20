import scvelo as scv
import joblib
import cellrank as cr
import pandas as pd

adata = scv.read('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/scvelo/scvelo_result.h5ad')

vk = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel.h5ad")

ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
g = cr.estimators.GPCCA(combined_kernel)

g.fit(cluster_key="class", n_states=[4, 12])
g.plot_macrostates(which="all", discrete=True, legend_loc="right", s=100, save='macrostates.png')

g.predict_terminal_states(allow_overlap=True)
g.plot_macrostates(which="terminal", legend_loc="right", s=100, save='macrostates_terminal.png')

g.plot_macrostates(which="terminal", discrete=False, save='macrostates_terminal_cont.png')

g.predict_initial_states(allow_overlap=True)
g.plot_macrostates(which="initial", legend_loc="right", s=100, save='macrostates_initial.png')

joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator.h5ad")
