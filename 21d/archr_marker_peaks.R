library(ArchR)
library(parallel)
#library(Seurat)
#library(Signac)
#library(BSgenome.Mmusculus.UCSC.mm10)

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all')
groupBy <- 'clusters0.1'
output <- '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/clusters_01_dar.csv'
#data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.QC.annotated.RDS')

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
# 
# project <- addUMAP(
#   ArchRProj = project, 
#   reducedDims = "IterativeLSI", 
#   name = "UMAP", 
#   nNeighbors = 30, 
#   minDist = 0.5, 
#   metric = "cosine"
# )
# 
# project <- addImputeWeights(project)

pathToMacs2 <- findMacs2()

# idxSample <- BiocGenerics::which(project$sex %in% "f")
# cellsSample <- project$cellNames[idxSample]
# proj.f <- project[cellsSample, ]
# 
# idxSample <- BiocGenerics::which(proj.f$clusters0.1 %in% "0")
# cellsSample <- proj.f$cellNames[idxSample]
# proj.f.c0 <- proj.f[cellsSample, ]

# Generic work flow for a project - repeated for subset projects above ---------
project <- addGroupCoverages(ArchRProj = project, groupBy = groupBy, force=TRUE)

project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = groupBy, 
  pathToMacs2 = pathToMacs2
)

project <- addPeakMatrix(project)

#saveArchRProject(project)

markersPeaks <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = groupBy,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.9 & Log2FC >= 0.1") # had to loosen cutoffs (default "FDR <= 0.01 $ Log2FC >= 1")
write.csv(markerList, output)
