import cellrank as cr
import scanpy as sc
import joblib

# https://cellrank.readthedocs.io/en/latest/notebooks/tutorials/estimators/800_gene_trends.html#visualize-expression-trends-via-line-plots
g = joblib.load("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/cellrank/gpcca_estimator_reint_stages.h5ad")

model = cr.models.GAM(g.adata)

genes = ['Cpne6', 'Celf4', 'Apc', 'Epha5', 'Atp1b1', 'Ptprn', 'Son', 'Atp2b1', 'AC149090.1', 'Gria4', 'Il1rapl1', 'Gad1', 'Ntm', 'Oxr1', 'Atp1a3', 'Snhg11', 'Srrm2', 'Eml5', 'Trdn', 'Btg1', 'Kmt2a', 'Cbfa2t3', 'Robo3', 'Itm2b', 'Gm32618', 'Nit1', 'Klhl26', 'Meis2', 'Plxna2', 'Pcp4', 'Ubb', 'Adam33', 'Bach2', 'Nrgn']

cr.pl.gene_trends(
    g.adata,
    model=model,
    data_key="velocity",
    genes=genes,
    same_plot=True,
    ncols=3,
    time_key="latent_time",
    hide_cells=True,
    save = "reint_clusters_gene_trends.png"
)
