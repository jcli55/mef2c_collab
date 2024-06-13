library(Seurat)
library(Signac)

# Script to recreate the Seurat object with just the RNA raw counts for linking with ArchR project (requires unmodified seurat object)
# Create the Seurat objects with raw counts and merge them

data.dir <- c('/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/f_con/filtered_feature_bc_matrix.h5', 
              '/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/f_ko/filtered_feature_bc_matrix.h5', 
              '/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/m_con/filtered_feature_bc_matrix.h5')

f_con_counts <- Read10X_h5(data.dir[1])
f_ko_counts <- Read10X_h5(data.dir[2])
m_con_counts <- Read10X_h5(data.dir[3])

f_con <- CreateSeuratObject(f_con_counts$'Gene Expression', project = 'f_con')
f_con <- RenameCells(f_con, new.names = paste('f_con', colnames(f_con), sep='_'))

f_ko <- CreateSeuratObject(f_ko_counts$'Gene Expression', project = 'f_ko')
f_ko <- RenameCells(f_ko, new.names = paste('f_ko', colnames(f_ko), sep='_'))

m_con <- CreateSeuratObject(m_con_counts$'Gene Expression', project = 'm_con')
m_con <- RenameCells(m_con, new.names = paste('m_con', colnames(m_con), sep='_'))

raw <- merge(f_con, list(f_ko, m_con))

# Filter out cells and transfer labels

data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.ATAC.QC.annotated.RDS')

intersect <- intersect(colnames(raw), colnames(data))

subset <- subset(raw, cells = intersect)

df <- data.frame(data@meta.data)
annotation <- c()

for (x in colnames(subset)) {
  annotation <- append(annotation, df[x, 'res0.5.annot'])
}

subset$res0.5.annot <- annotation

saveRDS(subset, '/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.QC.annotated.RDS')