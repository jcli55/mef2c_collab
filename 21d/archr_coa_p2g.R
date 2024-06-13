library(ArchR)
library(parallel)

# This script separates the data by genotype (con vs ko) and computes and plots co-accessibility and peak2gene linkage

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all_2/')

idxSample <- BiocGenerics::which(project$geno %in% "con")
cellsSample <- project$cellNames[idxSample]
con <- project[cellsSample, ]

idxSample <- BiocGenerics::which(project$geno %in% "ko")
cellsSample <- project$cellNames[idxSample]
ko <- project[cellsSample, ]

saveArchRProject(con, outputDirectory = "/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_con/", load = TRUE)
saveArchRProject(ko, outputDirectory = "/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_ko/", load = TRUE)

# First make sure peaks are called based on clusters0.1
con <- addGroupCoverages(ArchRProj = con, groupBy = "clusters0.1", force = TRUE)
ko <- addGroupCoverages(ArchRProj = ko, groupBy = "clusters0.1", force = TRUE)

pathToMacs2 <- findMacs2()

con <- addReproduciblePeakSet(
  ArchRProj = con, 
  groupBy = "clusters0.1", 
  pathToMacs2 = pathToMacs2
)
ko <- addReproduciblePeakSet(
  ArchRProj = ko, 
  groupBy = "clusters0.1", 
  pathToMacs2 = pathToMacs2
)

con <- addPeakMatrix(con)
ko <- addPeakMatrix(ko)

# Get marker genes
markerGenes  <- c('Eps8l1', 'Crhbp', 'Adam33', 'Nptx2', 'Cfap54', 'Cmya5', 'C1ql2', 'Sgcd', 'Sema3d', 'Dpy19l1')

# Co-accessibility
con <- addCoAccessibility(
  ArchRProj = con,
  reducedDims = "IterativeLSI"
)
ko <- addCoAccessibility(
  ArchRProj = ko,
  reducedDims = "IterativeLSI"
)

p <- plotBrowserTrack(
  ArchRProj = con, 
  groupBy = "clusters0.1", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(con, resolution = 1000)
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility-Con.pdf", 
        ArchRProj = con, 
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
  ArchRProj = ko, 
  groupBy = "clusters0.1", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(ko, resolution = 1000)
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility-Ko.pdf", 
        ArchRProj = ko, 
        addDOC = FALSE, width = 5, height = 5)

# Peaks-to-Genes
con <- addPeak2GeneLinks(
  ArchRProj = con,
  reducedDims = "IterativeLSI"
)
ko <- addPeak2GeneLinks(
  ArchRProj = ko,
  reducedDims = "IterativeLSI"
)

p <- plotBrowserTrack(
  ArchRProj = con, 
  groupBy = "clusters0.1", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(con, resolution = 1000)
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-Con.pdf", 
        ArchRProj = con, 
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
  ArchRProj = ko, 
  groupBy = "clusters0.1", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(ko, resolution = 1000)
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-Ko.pdf", 
        ArchRProj = ko, 
        addDOC = FALSE, width = 5, height = 5)

p2 <- plotPeak2GeneHeatmap(ArchRProj = con, groupBy = "clusters0.1")
p3 <- plotPeak2GeneHeatmap(ArchRProj = ko, groupBy = "clusters0.1")

plotPDF(plot = p2, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-Heatmap-Con.pdf", 
        ArchRProj = con, 
        addDOC = FALSE)
plotPDF(plot = p3, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-Heatmap-Ko.pdf", 
        ArchRProj = ko, 
        addDOC = FALSE)

saveArchRProject(con)
saveArchRProject(ko)