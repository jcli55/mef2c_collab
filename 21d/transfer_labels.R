library(ArchR)
library(parallel)
library(Seurat)
library(Signac)

# This script transfers metadata from the Seurat Object Suyang prepared (integration + annotation) to the ArchR project

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/merged/')
data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.ATAC.QC.annotated.RDS')

project$sampleid <- gsub('#','_', project$cellNames)
df <- DataFrame(data@meta.data)

intersect <- intersect(colnames(data), project$sampleid)

idxSample <- BiocGenerics::which(project$sampleid %in% intersect)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

clusters <- c()
annotation <- c()
sex <- c()
geno <- c()
clusters0.1 <- c()

for (x in project$sampleid) {
	clusters <- append(clusters, df[x, 'res0.5.reorder'])
	annotation <- append(annotation, df[x, 'res0.5.annot'])
	sex <- append(sex, df[x, 'sex'])
	geno <- append(geno, df[x, 'geno'])
	clusters0.1 <- append(clusters0.1, df[x, 'integrated_snn_res.0.1'])
}

project$clusters <- clusters
project$annotation <- annotation
project$sex <- sex
project$geno <- geno
project$clusters0.1 <- clusters0.1

saveArchRProject(project)