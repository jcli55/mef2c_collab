library(ArchR)
library(parallel)

# Plot p2g browser tracks for the new dataset

# Load the projects (created in archr_p2g.R)
project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/project')
projList <- c()
for(sam in unique(project$Sample)) {
  projList <- append(projList, loadArchRProject(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam)))
}

# Create a df of the genes and their lengths for the genes on the reverse strand
df <- data.frame(gene = c('Adam33','Crhbp','Dpy19l1','Cfap54','Pde10a','Bach2','Fosb','Bcl2','Eps8l1','Nptx2','Npas4','Itm2b','Mctp1','Htr2a','Numb','Slit2','Sgcd','Slc15a2','Sema3d','Tubb3','Matn2','Capn11','Cpne7'),
                 extend.up = c(50000, 50000, 140000, 350000, 100000, 100000, 10000, 200000, 50000, 10000, 10000, 40000, 300000, 50000, 150000, 100000, 1200000, 40000, 50000, 9000, 40000, 40000, 25000),
                 extend.down = c(50000, 50000, 50000, 50000, 600000, 400000, 15000, 120000, 50000, 15000, 15000, 20000, 750000, 90000, 20000, 450000, 200000, 20000, 250000, 25000, 170000, 40000, 50000),
                 regulation = c('down','down','down','down','down','down','down','down','down','down','down','down','down','down','down','up','up','up','up','up','up','up','up')
)

# Plot the browser tracks
for (proj in projList) {
  plotp2g.down <- c()
  plotp2g.up <- c()
  for (i in 1:length(df$gene)) {
    p <- plotBrowserTrack(
      ArchRProj = proj, 
      groupBy = "annot.res0.05", 
      geneSymbol = df$gene[i], 
      upstream = df$extend.up[i],
      downstream = df$extend.down[i],
      loops = getPeak2GeneLinks(proj, resolution = 1000)
    )
    
    if(df$regulation[i] == 'down') {
      plotp2g.down <- append(plotp2g.down, p)
    } else {
      plotp2g.up <- append(plotp2g.up, p)
    }
  }
  
  plotPDF(plotlist = plotp2g.down, 
          name = paste0("Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-", unique(proj$Sample), "-down.pdf"), 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
  plotPDF(plotlist = plotp2g.up, 
          name = paste0("Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-", unique(proj$Sample), "-up.pdf"), 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
}
