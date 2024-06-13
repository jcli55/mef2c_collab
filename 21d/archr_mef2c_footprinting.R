library(ArchR)
library(parallel)

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all')

# Subset the project
idxSample <- BiocGenerics::which(project$sex %in% "f")
cellsSample <- project$cellNames[idxSample]
proj.f <- project[cellsSample, ]

idxSample <- BiocGenerics::which(proj.f$geno %in% "con")
cellsSample <- proj.f$cellNames[idxSample]
proj.f.con <- proj.f[cellsSample, ]

idxSample <- BiocGenerics::which(proj.f$geno %in% "ko")
cellsSample <- proj.f$cellNames[idxSample]
proj.f.ko <- proj.f[cellsSample, ]

# Project f_con
for (proj in c(proj.f.con, proj.f.ko)) {
  motifPositions <- getPositions(proj)
  
  motifs <- c("Mef2c")
  markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
  markerMotifs
  
  proj <- addGroupCoverages(ArchRProj = proj, groupBy = "clusters0.1", force = TRUE)
  
  seFoot <- getFootprints(
    ArchRProj = proj, 
    positions = motifPositions[markerMotifs], 
    groupBy = "clusters0.1"
  )
  
  plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj, 
    normMethod = "Subtract",
    plotName = paste0("Footprints-F-", unique(proj$geno)),
    addDOC = FALSE,
    smoothWindow = 5
  )
}