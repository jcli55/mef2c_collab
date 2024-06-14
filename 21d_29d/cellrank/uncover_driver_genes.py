import cellrank as cr
import scanpy as sc
import joblib

#adata = sc.read_h5ad('/storage/singlecell/jeanl/organoid/data/Chen/h5ad/multivelo/multivelo_result.h5ad') # Will need if plan to plot some embeddings

g = joblib.load("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/gpcca_estimator_multivelo.h5ad")

driver_df = g.compute_lineage_drivers()
driver_df.to_csv("/storage/singlecell/jeanl/organoid/csv/lineage_drivers.csv")

# More plots and functionality in the tutorial