library(ArchR)
library(parallel)

# Motif deviations (13.1)
# Note: if errors occur during addBgdPeaks() or addDeviationsMatrix(), then try deleting the GroupCoverages and PeakCalls directories in the saved project and rerun GroupCoverages and Peak Calling

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all')

# Make sure Motif annotation was run
if("Motif" %ni% names(project@peakAnnotation)){
  project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif")
}

# Add Motif Deviations
project <- addBgdPeaks(project)

project <- addDeviationsMatrix(
  ArchRProj = project, 
  peakAnnotation = "Motif",
  force = TRUE
)

saveArchRProject(project)

# Can plot the motif deviation scores
#plotVarDev <- getVarDeviations(project, name = "MotifMatrix", plot = TRUE)
#plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = project, addDOC = FALSE)

# Plot motif z-scores on UMAP
motifs <- c('Mef2c','Mef2a','Meis1','Meis2','Bach2','Tcf12','Tcf4','Klf12','Sox5','Pbx1','Pbx3','Runx2','Pparg','Rora','Nfib','Nfia','Fosb','Junb')
markerMotifs <- getFeatures(project, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(project)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# Plot gene score on UMAP
markerRNA <- getFeatures(project, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- markerRNA[markerRNA %ni% c('Bach2os', 'Ppargc1a', 'Ppargc1b', 'Sox5os3')]

p3 <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(project)
)

p4 <- lapply(p3, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p4))

# Plot gene integration on UMAP
markerRNA <- getFeatures(project, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
markerRNA <- markerRNA[markerRNA %ni% c('Bach2os', 'Ppargc1a', 'Ppargc1b', 'Sox5os3')]

p5 <- plotEmbedding(
  ArchRProj = project, 
  colorBy = "GeneIntegrationMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP",
  continuousSet = "blueYellow",
  imputeWeights = getImputeWeights(project)
)

p6 <- lapply(p5, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p6))

plotPDF(p2,p4,p6, name = "Plot-UMAP-Motif-Deviations.pdf", ArchRProj = project, addDOC = FALSE, width = 5, height = 5)
