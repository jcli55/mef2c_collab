import scanpy as sc

path = "/storage/chentemp/suyangb/m2c_multiome/data/"
samples = ["M2C_21d_f_con",
            "M2C_21d_f_ko",
            "M2C_21d_m_con",
            "M2C_21d_m_ko_repeat_final",
            "M2C_29d_m_con_final",
            "M2C_29d_m_ko_final"]
id = ["21d_f_con",
        "21d_f_ko",
        "21d_m_con",
        "21d_m_ko",
        "29d_m_con",
        "29d_m_ko"]
list = []

for i in range(0, len(samples)):
    temp = sc.read_10x_mtx(f'{path}{samples[i]}/filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True, gex_only=True)

    names = []
    for name in temp.obs_names:
        names.append(f'{id[i]}#{name}')
    temp.obs_names = names
    list.append(temp)

adata_rna = sc.concat(list)

adata_rna.write_h5ad("/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/multivelo/adata_rna.h5ad")
