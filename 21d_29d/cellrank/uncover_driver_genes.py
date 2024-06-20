import cellrank as cr
import scanpy as sc
import joblib

g = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator.h5ad")

driver_df = g.compute_lineage_drivers()
driver_df.to_csv("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/lineage_drivers.csv")

# More plots and functionality in the tutorial
