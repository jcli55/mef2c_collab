library(Seurat)
library(Signac)

# QC based on RNA assay
# Read in the raw seurat objects
data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.ATAC.RDS')

# Plot QC metrics and filter out low-quality cells
for (i in data) {
    VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    i <- subset(i, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 1.5))
}

# Save the filtered seurat objects
saveRDS(data, '/storage/chentemp/u250758/mef2c_collab/data/rds/qc/seurat.objs.RNA.QC.RDS')

# ----- QC based on ATAC assay using ArchR (see archr_ATAC_qc.R) -----
# Load the cells that passed QC from ArchR (reformat to match seurat cellnames), and keep only intersecting cells
for (i in data) {
    archr.cells <- read.csv(paste0('/storage/chentemp/u250758/mef2c_collab/data/archr/csv/', i, '_ATAC_qc_cells.csv'))
    archr.cells$x <- gsub('#', '_', archr_cells$x)
    i <- subset(i, cells = intersect(colnames(i), archr.cells$x))
}

saveRDS(data, '/storage/chentemp/u250758/mef2c_collab/data/rds/qc/seurat.objs.RNA.ATAC.QC.RDS')