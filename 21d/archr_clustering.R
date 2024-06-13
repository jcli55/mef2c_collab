library(ArchR)
library(parallel)

# This script adds clusters based on the dimensional reduction 'IterativeLSI' based on the ATAC data

# Load the merged ArchR project (created with 'create_archr_project.R' and then 'transfer_labels.R')
project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/merged')

# addIterativeLSI() and addUMAP() already ran in 'archr_marker_peaks.R')
# project <- addIterativeLSI(
#   ArchRProj = project,
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI", 
#   iterations = 2, 
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.2), 
#     sampleCells = 10000, 
#     n.start = 10
#   ), 
#   varFeatures = 25000, 
#   dimsToUse = 1:30
# )

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ArchR.Clusters.0.8",
  resolution = 0.8 #Change the resolution here
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ArchR.Clusters.1.0",
  resolution = 1.0
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ArchR.Clusters.0.5",
  resolution = 0.5
)

project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ArchR.Clusters.0.1",
  resolution = 0.1
)

# project <- addUMAP(
#   ArchRProj = project, 
#   reducedDims = "IterativeLSI", 
#   name = "UMAP", 
#   nNeighbors = 30, 
#   minDist = 0.5, 
#   metric = "cosine"
# )

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = c("Sample", "annotation", "geno", "clusters", "clusters0.1"), embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ArchR.Clusters.1.0", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ArchR.Clusters.0.8", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ArchR.Clusters.0.5", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "ArchR.Clusters.0.1", embedding = "UMAP")
plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(project)