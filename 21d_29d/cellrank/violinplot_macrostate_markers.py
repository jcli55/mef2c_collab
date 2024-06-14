import scanpy as sc

adata = sc.read_h5ad("/storage/singlecell/jeanl/organoid/data/Chen/h5ad/cellrank/velocity_kernel_multivelo.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# compute a score in scanpy by aggregating across a few ductal markers
sc.tl.score_genes(adata, gene_list=['SFRP2','NFIA','ASCL1'], score_name="initial_score", use_raw=False)

sc.tl.score_genes(adata, gene_list=['TFAP2A','GAD1','GAD2'], score_name="AC_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['VSX1','VSX2'], score_name="BC_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['ONECUT1','ONECUT2'], score_name="HC_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['NEFM','RBPMS','POU4F2'], score_name="RGC_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['RLBP1','SLC1A3','SFRP2'], score_name="MG_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['ARR3'], score_name="Cone_score", use_raw=False)
sc.tl.score_genes(adata, gene_list=['RCVRN','CRX','PROM1'], score_name="Rod_score", use_raw=False) # No RHO

# visualize via heatmaps
sc.pl.violin(adata, keys="initial_score", groupby="macrostates", rotation=90, save='macrostates_prpc_markers.png', use_raw=False)
sc.pl.violin(adata, keys="AC_score", groupby="macrostates", rotation=90, save='macrostates_AC_markers.png', use_raw=False)
sc.pl.violin(adata, keys="BC_score", groupby="macrostates", rotation=90, save='macrostates_BC_markers.png', use_raw=False)
sc.pl.violin(adata, keys="HC_score", groupby="macrostates", rotation=90, save='macrostates_HC_markers.png', use_raw=False)
sc.pl.violin(adata, keys="RGC_score", groupby="macrostates", rotation=90, save='macrostates_RGC_markers.png', use_raw=False)
sc.pl.violin(adata, keys="MG_score", groupby="macrostates", rotation=90, save='macrostates_MG_markers.png', use_raw=False)
sc.pl.violin(adata, keys="Cone_score", groupby="macrostates", rotation=90, save='macrostates_Cone_markers.png', use_raw=False)
sc.pl.violin(adata, keys="Rod_score", groupby="macrostates", rotation=90, save='macrostates_Rod_markers.png', use_raw=False)