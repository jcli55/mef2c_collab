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

# Load in the DEGs
degDo <- c('Eps8l1','Crhbp','Adam33','Nptx2','Cfap54','Cmya5','Bub1b','Dpy19l1','Mid1','Rasgef1c','Csmd2','Rgs6','Tpbg')
degUp <- c('C1ql2','Capn11','Sgcd','Sema3d','Rec8','Wnt9b','Slit2','Th','Syt5','Gapdh','Ndufa4','Tuba1a','Itm2c','Ubb')

# Plot the gene scores
for (proj in c(proj.f.con, proj.f.ko)) {
  geno <- unique(proj$geno)
  count <- 0
  proj <- addImputeWeights(proj)
  for (geneList in list(degDo, degUp)) {
    count <- count + 1
    
    p <- plotEmbedding(
      ArchRProj = proj, 
      colorBy = "GeneScoreMatrix", 
      name = geneList, 
      embedding = "UMAP",
      imputeWeights = getImputeWeights(proj)
    )
    
    p2 <- lapply(p, function(x){
      x + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 3) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()
        )
    })
    #do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
    
    plotPDF(do.call(cowplot::plot_grid, c(list(ncol = 3),p2)), 
            name = paste0("Plot-UMAP-", count, "-Genes-W-Imputation-", geno, ".pdf"), 
            ArchRProj = proj, 
            addDOC = FALSE, width = 5, height = 5)
  }
}

# Get the average gene scores for each gene between f_con_all vs f_ko_all and f_con_c0 vs f_ko_c0
idxSample <- BiocGenerics::which(proj.f$clusters0.1 %in% "0")
cellsSample <- proj.f$cellNames[idxSample]
proj.f.c0 <- proj.f[cellsSample, ]

seGS.all <- getGroupSE(proj.f, useMatrix = "GeneScoreMatrix", groupBy = "geno")
seGS.c0 <- getGroupSE(proj.f.c0, useMatrix = "GeneScoreMatrix", groupBy = "geno")

rownames(seGS.all) <- rowData(seGS.all)$name
rownames(seGS.c0) <- rowData(seGS.c0)$name

write.csv(assays(seGS.all)$GeneScoreMatrix, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/avg_genescore_all.csv')
write.csv(assays(seGS.c0)$GeneScoreMatrix, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/avg_genescore_c0.csv')