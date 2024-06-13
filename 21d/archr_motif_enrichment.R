library(ArchR)
library(parallel)

# Motif enrichment

proj.f <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_all')
prof.f.c0 <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/proj_f_c0')

# Already ran in "archr_marker_peaks.R"
markersPeaks.f <- getMarkerFeatures(
  ArchRProj = proj.f, 
  useMatrix = "PeakMatrix", 
  groupBy = "geno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks.f.c0 <- getMarkerFeatures(
  ArchRProj = proj.f.c0, 
  useMatrix = "PeakMatrix", 
  groupBy = "geno",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

proj.f <- addMotifAnnotations(ArchRProj = prof.f, motifSet = "cisbp", name = "Motif")
proj.f.c0 <- addMotifAnnotations(ArchRProj = prof.f.c0, motifSet = "cisbp", name = "Motif")

# Enriched in f con all --------------------------------------------------------
motifsUp.f <- peakAnnoEnrichment(
  seMarker = markersPeaks.f,
  ArchRProj = proj.f,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp.f), mlog10Padj = assay(motifsUp.f)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.csv(df, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/f_all_motifs_up.csv')

ggUp.all <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggUp.all, name = "Plot-Motifs-Enriched-con-all.pdf", ArchRProj = proj.f, addDOC = FALSE, width = 5, height = 5)

# Enriched in f con GC ---------------------------------------------------------
motifsUp.f.c0 <- peakAnnoEnrichment(
  seMarker = markersPeaks.f.c0,
  ArchRProj = proj.f.c0,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp.f.c0), mlog10Padj = assay(motifsUp.f.c0)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.csv(df, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/f_c0_motifs_up.csv')

ggUp.c0 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggUp.c0, name = "Plot-Motifs-Enriched-con-c0.pdf", ArchRProj = proj.f, addDOC = FALSE, width = 5, height = 5)

# Enriched in f ko all ---------------------------------------------------------
motifsDo.f <- peakAnnoEnrichment(
  seMarker = markersPeaks.f,
  ArchRProj = proj.f,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsDo.f), mlog10Padj = assay(motifsDo.f)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.csv(df, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/f_all_motifs_do.csv')

ggDo.all <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggDo.all, name = "Plot-Motifs-Enriched-ko-all.pdf", ArchRProj = proj.f, addDOC = FALSE, width = 5, height = 5)

# Enriched in f ko GC ----------------------------------------------------------
motifsDo.f.c0 <- peakAnnoEnrichment(
  seMarker = markersPeaks.f.c0,
  ArchRProj = proj.f.c0,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 1 & Log2FC <= -0" # loosened criteria
)

df <- data.frame(TF = rownames(motifsDo.f.c0), mlog10Padj = assay(motifsDo.f.c0)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.csv(df, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/f_c0_motifs_do.csv')

ggDo.c0 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black",
    max.overlaps = 100
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggDo.c0, name = "Plot-Motifs-Enriched-ko-c0.pdf", ArchRProj = proj.f, addDOC = FALSE, width = 5, height = 5)