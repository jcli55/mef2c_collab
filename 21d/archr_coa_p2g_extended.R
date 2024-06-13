library(ArchR)
library(parallel)

# This script extends the browser tracks window to cover large gene bodies (Cfap54, Cmya5, Sgcd, Sema3d, Dpy19l1)
# Extension from archr_coa_p2g.R

con <- loadArchRProject("/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_con/")
ko <- loadArchRProject("/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_ko/")

# Create a df of the genes and their lengths for the genes on the reverse strand
df <- data.frame(gene = c('Cfap54', 'Cmya5', 'Dpy19l1'),
		            length = c(305418, 104012, 91363))

# Plot the browser tracks
for (project in c(con, ko)) {
  plotcoa <- c()
  plotp2g <- c()
	for (i in 1:length(df$gene)) {
	  p <- plotBrowserTrack(
	    ArchRProj = project, 
	    groupBy = "clusters0.1", 
	    geneSymbol = df$gene[i], 
	    upstream = df$length[i]+50000,
	    downstream = 50000,
	    loops = getCoAccessibility(project, resolution = 1000)
	  )
	  plotcoa <- append(plotcoa, p)
	  
	  p2 <- plotBrowserTrack(
	    ArchRProj = project, 
	    groupBy = "clusters0.1", 
	    geneSymbol = df$gene[i], 
	    upstream = df$length[i]+50000,
	    downstream = 50000,
	    loops = getPeak2GeneLinks(project, resolution = 1000)
	  )
	  plotp2g <- append(plotp2g, p2)
	}

  # Wanted to extend the browser window for Sgcd 200,000bp instead of 50,000
  p <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "clusters0.1", 
    geneSymbol = 'Sgcd', 
    upstream = 1000000+200000,
    downstream = 200000,
    loops = getCoAccessibility(project, resolution = 1000)
  )
  plotcoa <- append(plotcoa, p)
  
  p2 <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "clusters0.1", 
    geneSymbol = 'Sgcd', 
    upstream = 1000000+200000,
    downstream = 200000,
    loops = getPeak2GeneLinks(project, resolution = 1000)
  )
  plotp2g <- append(plotp2g, p2)

  # Now add Sema3d (which is not on the reverse strand and needs to extend downstream)
  p <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "clusters0.1", 
    geneSymbol = 'Sema3d', 
    upstream = 50000,
    downstream = 205564+50000,
    loops = getCoAccessibility(project, resolution = 1000)
  )
  plotcoa <- append(plotcoa, p)
  
  p2 <- plotBrowserTrack(
    ArchRProj = project, 
    groupBy = "clusters0.1", 
    geneSymbol = 'Sema3d', 
    upstream = 50000,
    downstream = 205564+50000,
    loops = getPeak2GeneLinks(project, resolution = 1000)
  )
  plotp2g <- append(plotp2g, p2)
  
  plotPDF(plotlist = plotcoa, 
          name = paste0("Plot-Tracks-Marker-Genes-with-CoAccessibility-", unique(project$geno), "-extended.pdf"), 
          ArchRProj = project, 
          addDOC = FALSE, width = 5, height = 5)
  plotPDF(plotlist = plotp2g, 
          name = paste0("Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-", unique(project$geno), "-extended.pdf"), 
          ArchRProj = project, 
          addDOC = FALSE, width = 5, height = 5)
}