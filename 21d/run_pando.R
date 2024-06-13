library(tidyverse)
library(Pando)
library(Seurat)
library(Signac)
library(doParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)

data('motifs')
data('motif2tf')

# This script runs the basic Pando workflow
# https://quadbio.github.io/Pando/articles/getting_started.html

# For loop to run the workflow through each sample - samples created with 'prepare_pando_objects.RDS'
samples <- c('f.con', 'f.con.c0', 'f.ko', 'f.ko.c0')

# Function to capitalize the first letter of a string while lower-casing every other letter
# Used to match formatting between mm10 genome and Pando's curated motif list
# Source: https://github.com/quadbio/Pando/issues/30#issue-1564838244 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
  x
}

motif2tf$tf <- firstup(motif2tf$tf)

for (sam in samples) {
  # Load in the data and normalize it (Normalizing steps from these vignettes:)
  # https://satijalab.org/seurat/articles/pbmc3k_tutorial
  # https://stuartlab.org/signac/articles/pbmc_vignette 
  data <- c()
  data <- readRDS(paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.', sam, '.RDS'))
  
  # Peaks on 'GL456216.1','GL456233.1','JH584292.1','JH584304.1','JH584295.1' and reference mm10 genome chrM are clashing
  # Removing chrM from annotation:
  
  annotation <- data@assays$ATAC@annotation
  annotation <- annotation[annotation@seqnames != 'chrM']
  seqlevels(annotation) <- seqlevels(annotation)[seqlevels(annotation) != 'chrM']
  data@assays$ATAC@annotation <- annotation
  
  DefaultAssay(data) <- 'RNA'
  data <- NormalizeData(data) 
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(data) 
  data <- ScaleData(data, features = all.genes)
  
  DefaultAssay(data) <- 'ATAC'
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 'q0')
  data <- RunSVD(data)
  
  # Start the Pando workflow ---------------------------------------------------
  # Initiate GRN
  data <- initiate_grn(
    data,
    rna_assay = 'RNA',
    peak_assay = 'ATAC'
  )
  
  # Scanning for motifs
  data <- find_motifs(
    data, 
    pfm = motifs, 
    motif_tfs = motif2tf,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )
  
  # Inferring the GRN
  registerDoParallel(4)
  data <- infer_grn(
    data,
    peak_to_gene_method = 'GREAT',
    genes = data@assays$RNA@var.features,
    parallel = T,
    tf_cor = 0.001,
    upstream = 1e+05,
    downstream = 1e+05,
    extend = 1e+06
  )
  
  saveRDS(data, paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/', sam, '.grn.RDS'))
  write.csv(coef(data), paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/', sam, '_infer_grn.csv'))
  
  # Module discovery
  # data <- find_modules(data) #default params
#   data <- find_modules( #these are more lenient params
#     data,
#     p_thresh = 1,
#     nvar_thresh = 1,
#     min_genes_per_module = 1,
#     rsq_thresh = 0
#   )
#   #  data <- find_modules( #these are even more lenient params I had to use for Rod ro
#   #    data, 
#   #    p_thresh = 0.2,
#   #    nvar_thresh = 2, 
#   #    min_genes_per_module = 1, 
#   #    rsq_thresh = 0.025
#   #  )
#   
#   # Some other plots from the vignette
#   #plot_gof(data, point_size=3)
#   #plot_module_metrics(data)
#   
#   # Visualizing the GRN (The png + dev.off doesn't work in scripts)
#   #data <- get_network_graph(data)
#   #png(filename = paste0('/storage/singlecell/jeanl/organoid/figures/grns/', mclass, '_grn.png'))
#   #plot_network_graph(data)
#   #dev.off()
#   
#   write.csv(NetworkModules(data)@meta, paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/', sam, '.grn.csv'))
}
