library(ArchR)
library(parallel)
addArchRThreads(threads = 8)
addArchRGenome('mm10')

# Setup the input files, sample names, and output files
inputFiles <- c('/storage/chentemp/u250758/mef2c_collab/data/29d/cellranger_outs/m_con/outs/atac_fragments.tsv.gz',
				'/storage/chentemp/u250758/mef2c_collab/data/29d/cellranger_outs/m_ko/outs/atac_fragments.tsv.gz')
names <- c( "m_con","m_ko") 

#please modify accordingly
parent_dir <- '.../29d/....'

outputFiles <- c(
		 		'/storage/chentemp/u250758/mef2c_collab/data/29d/archr/projects/proj_m_con',
		 		'/storage/chentemp/u250758/mef2c_collab/data/29d/archr/projects/proj_m_ko')

outputFilesQC <- c(
				'/storage/chentemp/u250758/mef2c_collab/data/29d/archr/projects/filter_doublets/proj_m_con',
		 		'/storage/chentemp/u250758/mef2c_collab/data/29d/archr/projects/filter_doublets/proj_m_ko')

#----
# filterFrags is no longer a valid input. Please use minFrags! Setting filterFrags value to minFrags!
#   filterTSS is no longer a valid input. Please use minTSS! Setting filterTSS value to minTSS!
#   Error in names[i] : object of type 'builtin' is not subsettable

# For each sample----
for (i in 1:2) {
	# Create the ArchR project from the fragments.tsv.gz files
	ArrowFiles <- createArrowFiles(
	  inputFiles = inputFiles[i],
	  sampleNames = names[i],
	  filterTSS = 4, #Dont set this too high because you can always increase later
	  filterFrags = 1000, 
	  addTileMat = TRUE,
	  addGeneScoreMat = TRUE
	)

	doubScores <- addDoubletScores(
	  input = ArrowFiles,
	  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
	  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
	  LSIMethod = 1
	)

	proj <- ArchRProject(
	  ArrowFiles = ArrowFiles, 
	  outputDirectory = outputFiles[i],
	  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
	)
	
	# Save the raw project
	saveArchRProject(ArchRProj = proj, outputDirectory = outputFiles[i], load = FALSE)
	
	# Filter out doublets and save the filtered project and cells that passed QC
	proj <- filterDoublets(proj)

	saveArchRProject(ArchRProj = proj, outputDirectory = outputFilesQC[i], load = FALSE)
	write.csv(proj$cellNames, paste0('/storage/chentemp/u250758/mef2c_collab/data/29d/archr/csv/', names[i], '_ATAC_qc_cells.csv'))
}