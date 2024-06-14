import cellrank as cr
import scanpy as sc
import joblib

g = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo.h5ad")
g = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo_scvi_states.h5ad")

g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save='fate_probabilities_sep.svg')
g.plot_fate_probabilities(same_plot=True, save='fate_probabilities.png')

joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo.h5ad")
joblib.dump(g, "/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo_scvi_states.h5ad")