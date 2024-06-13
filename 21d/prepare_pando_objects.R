library(Seurat)
library(Signac)

# Prepare the objects for Pando pipeline
data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.QC.annotated.RDS') # Created with 'create_merge_seurat_obj.R'
merged.data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.ATAC.QC.annotated.RDS') # Object created by Suyang

metadata <- merged.data@meta.data # Also saved in '/storage/chentemp/u250758/mef2c_collab/data/rds/metadata.csv'
intersect <- intersect(metadata$X, colnames(data))
data <- subset(data, cells = intersect)

data[['ATAC']] <- merged.data[['ATAC']]

saveRDS(data, '/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.ATAC.for.pando.RDS')

# Now prepare the subset objects for each sample:
# Get the female cells only
subset <- subset(data, cells = subset(metadata, sex =='f')$X) 

# Get the f_con cells only
subset <- subset(subset, cells = subset(metadata, geno == 'con')$X) 
saveRDS(subset, '/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.f.con.RDS')

# Get the f_ko cells only
subset <- subset(subset, cells = subset(metadata, geno == 'ko')$X) 
saveRDS(subset, '/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.f.ko.RDS')

# Get the f_con_c0 cells only
subset <- subset(subset, cells = subset(metadata, geno == 'con' & integrated_snn_res.0.1 == '0')$X) 
saveRDS(subset, '/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.f.con.c0.RDS')

# Get the f_ko_c0 cells only
subset <- subset(subset, cells = subset(metadata, geno == 'ko' & integrated_snn_res.0.1 == '0')$X) 
saveRDS(subset, '/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.f.ko.c0.RDS')