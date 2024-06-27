import cellrank as cr
import scanpy as sc
import joblib

# https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/estimators/800_gene_trends.html#visualize-expression-trends-via-line-plots
g = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator_orig_stage.h5ad")

model = cr.models.GAM(g.adata)

genes = ['Cpne6','Snhg11','Celf4','Apc','Oxr1','Ptprn','Son','Epha5','AC149090.1','Kcnh3','Trdn','Lrp8','Tcf4','Igsf9b','Taf6l','Adam33','Robo3','Bach2','Nrgn']

cr.pl.gene_trends(
    g.adata,
    model=model,
    data_key="velocity",
    genes=genes,
    same_plot=True,
    ncols=3,
    time_key="latent_time",
    hide_cells=True,
    save = "orig_clusters_gene_trends.png"
)
