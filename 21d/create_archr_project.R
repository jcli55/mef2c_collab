library(ArchR)
library(parallel)
addArchRThreads(threads = 8)
addArchRGenome('mm10')

# Creates the merged ArchR project for downstream analysis

# Setup the input files, sample names, and output files
inputFiles <- c('/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/f_con/atac_fragments.tsv.gz',
				'/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/f_ko/atac_fragments.tsv.gz',
				'/storage/chentemp/u250758/mef2c_collab/data/cellranger_outs/m_con/atac_fragments.tsv.gz')

names <- c('f_con', 'f_ko', 'm_con')

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names,
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#  LSIMethod = 1
#)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/storage/chentemp/u250758/mef2c_collab/data/archr/projects/merged',
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = '/storage/chentemp/u250758/mef2c_collab/data/archr/projects/merged', load = FALSE)
