library(ArchR)
library(parallel)

# Motif enrichment

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all')
project$clusters0.1 <- as.integer(project$clusters0.1)

# First make sure peaks are called based on the desired cellColData (clusters0.1 in this case)
pathToMacs2 <- findMacs2()

project <- addGroupCoverages(ArchRProj = project, groupBy = "clusters0.1", force = TRUE)

project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = "clusters0.1", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)
project <- addPeakMatrix(project, force = TRUE)

project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif", force = TRUE)

# Get the DARs
dar1 <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "clusters0.1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get the motif enrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker = dar1,
  ArchRProj = project,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-01", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)

# ------------------------------------------------------------------------------
# Same process but for clusters0.2

project <- addGroupCoverages(ArchRProj = project, groupBy = "clusters0.2", force = TRUE)

project <- addReproduciblePeakSet(
  ArchRProj = project, 
  groupBy = "clusters0.2", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)
project <- addPeakMatrix(project, force = TRUE)

project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif", force = TRUE)

dar2 <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "clusters0.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs2 <- peakAnnoEnrichment(
  seMarker = dar2,
  ArchRProj = project,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs2, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-02", width = 8, height = 6, ArchRProj = project, addDOC = FALSE)
