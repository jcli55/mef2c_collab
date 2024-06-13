library(ArchR)
library(parallel)
library(Seurat)

# Integrating RNA seq data to ATAC and linking peaks to genes

project <- loadArchRProject('/storage/chentemp/u250758/mef2c_collab/data/archr/projects/merged/')
data <- readRDS('/storage/chentemp/u250758/mef2c_collab/data/rds/seurat.objs.RNA.QC.annotated.RDS')

# Assume already ran archr_marker_peaks.R

groupList <- SimpleList(
  GC4.Lama3.5T4 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC4.Lama3.5T4')],
    RNA = colnames(data)[grep('GC4.Lama3.5T4', data$res0.5.annot)]
  ),
  GC8.Eda = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC8.Eda')],
    RNA = colnames(data)[grep('GC8.Eda', data$res0.5.annot)]
  ),
  GC9.active = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC9.active')],
    RNA = colnames(data)[grep('GC9.active', data$res0.5.annot)]
  ),
  PGC20.Vwc2l.Calb1 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('PGC20.Vwc2l.Calb1')],
    RNA = colnames(data)[grep('PGC20.Vwc2l.Calb1', data$res0.5.annot)]
  ),
  GC1 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC1')],
    RNA = colnames(data)[grep('GC1$', data$res0.5.annot)]
  ),
  PGC11.Kcnmb2.Calb2 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('PGC11.Kcnmb2.Calb2')],
    RNA = colnames(data)[grep('PGC11.Kcnmb2.Calb2', data$res0.5.annot)]
  ),
  C17.Grm3 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C17.Grm3')],
    RNA = colnames(data)[grep('C17.Grm3', data$res0.5.annot)]
  ),
  PGC16.Stac.Th = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('PGC16.Stac.Th')],
    RNA = colnames(data)[grep('PGC16.Stac.Th', data$res0.5.annot)]
  ),
  GC2 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC2')],
    RNA = colnames(data)[grep('GC2$', data$res0.5.annot)]
  ),
  C18.Grm3 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C18.Grm3')],
    RNA = colnames(data)[grep('C18.Grm3', data$res0.5.annot)]
  ),
  GC3 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC3')],
    RNA = colnames(data)[grep('GC3$', data$res0.5.annot)]
  ),
  PGC15.Kcnmb2.Calb2 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('PGC15.Kcnmb2.Calb2')],
    RNA = colnames(data)[grep('PGC15.Kcnmb2.Calb2', data$res0.5.annot)]
  ),
  GC6.Calb2 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC6.Calb2')],
    RNA = colnames(data)[grep('GC6.Calb2', data$res0.5.annot)]
  ),
  C14.Ndnf = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C14.Ndnf')],
    RNA = colnames(data)[grep('C14.Ndnf', data$res0.5.annot)]
  ),
  GC5.Kcnh8 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC5.Kcnh8')],
    RNA = colnames(data)[grep('GC5.Kcnh8', data$res0.5.annot)]
  ),
  GC7.Ak5= SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC7.Ak5')],
    RNA = colnames(data)[grep('GC7.Ak5', data$res0.5.annot)]
  ),
  C10.Ndnf = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C10.Ndnf')],
    RNA = colnames(data)[grep('C10.Ndnf', data$res0.5.annot)]
  ),
  C19.Cdh23 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C19.Cdh23')],
    RNA = colnames(data)[grep('C19.Cdh23', data$res0.5.annot)]
  ),
  GC13.transit.Kcnt2 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC13.transit.Kcnt2')],
    RNA = colnames(data)[grep('GC13.transit.Kcnt2', data$res0.5.annot)]
  ),
  C21 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C21')],
    RNA = colnames(data)[grep('C21$', data$res0.5.annot)]
  ),
  C22 = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('C22')],
    RNA = colnames(data)[grep('C22$', data$res0.5.annot)]
  ),
  GC12.transit = SimpleList(
    ATAC = project$cellNames[project$annotation %in% c('GC12.transit')],
    RNA = colnames(data)[grep('GC12.transit', data$res0.5.annot)]
  ),
  C23.24.25 = SimpleList( # Merged due to low cell counts
    ATAC = project$cellNames[project$annotation %in% c('C23.minor', 'C24', 'C25.Vip')],
    RNA = colnames(data)[grep('C23.minor|C24|C25.Vip', data$res0.5.annot)]
  )
)

project <- addGeneIntegrationMatrix(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = data,
  addToArrow = TRUE, # Set this to TRUE when satisfied with the labeling
  groupList = groupList,
  groupRNA = "res0.5.annot",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

project <- addPeak2GeneLinks(
  ArchRProj = project,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = project,
  corCutOff = 0.45,
  resolution = 1000, # Better for Browser Tracks
  returnLoops = TRUE
)

write.csv(p2g, '/storage/chentemp/u250758/mef2c_collab/data/archr/csv/p2g_GR.csv')

saveArchRProject(project)