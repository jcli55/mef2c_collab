library(ArchR)
library(parallel)

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/project')

# Make sure IterativeLSI is run
#project <- addIterativeLSI(
#  ArchRProj = project,
#  useMatrix = "TileMatrix", 
#  name = "IterativeLSI", 
#  iterations = 2, 
#  clusterParams = list( #See Seurat::FindClusters
#    resolution = c(0.2), 
#    sampleCells = 10000, 
#    n.start = 10
#  ), 
#  varFeatures = 25000, 
#  dimsToUse = 1:30
#)
#saveArchRProject('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/project')

# Subset based on sample
projList <- c()
for(sam in unique(project$Sample)) {
  idxSample <- BiocGenerics::which(project$Sample %in% sam)
  cellsSample <- project$cellNames[idxSample]
  proj <- project[cellsSample, ]
  
  saveArchRProject(proj, paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam))
  
  projList <- append(projList, loadArchRProject(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam)))
}

# Peak calling
pathToMacs2 <- findMacs2()
temp <- c()
for(proj in projList) {
  proj <- addGroupCoverages(ArchRProj = proj, groupBy = "annot.res0.05")
  
  proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "annot.res0.05", 
    pathToMacs2 = pathToMacs2
  )
  
  proj <- addPeakMatrix(proj)
  
  saveArchRProject(proj)
  
  temp <- append(temp, loadArchRProject(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam)))
}
projList <- temp

# P2G links and plot heatmap
plotList <- c()
for(proj in projList) {
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix"
  )
  
  p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "annot.res0.05")
  plotList <- append(plotList, p)
  
  saveArchRProject(proj)
}
plotPDF(plotList = plotList, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-Heatmap.pdf", 
        ArchRProj = project, 
        addDOC = FALSE)