import cellrank as cr
import scanpy as sc
import joblib

g = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator.h5ad")

g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False, save='fate_probabilities_sep.png')
g.plot_fate_probabilities(same_plot=True, save='fate_probabilities.png')

joblib.dump(g, "/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator.h5ad")
