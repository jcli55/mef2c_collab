library(tidyverse)
library(Pando)
library(Seurat)
library(Signac)
library(doParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)

# First create the mouse motif2tf table
# Load the data into R
url <- "https://jaspar2020.genereg.net/download/data/2020/CORE/JASPAR2020_CORE_non-redundant_pfms_jaspar.txt"
input_data <- readLines(url)

# Display the first few rows of the data
head(input_data)

# Read input data as a text connection
con <- textConnection(input_data)

# Create an empty data frame
df <- data.frame()

# Read lines from the text connection
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  # Check if the line starts with ">"
  if (substr(line, 1, 1) == ">") {
    # Extract the first and second columns from the line
    columns <- strsplit(line, "\\s+")[[1]]
    
    # Remove ">" from the first column
    columns[1] <- substr(columns[1], 2, nchar(columns[1]))
    
    # Append the columns to the data frame
    df <- rbind(df, data.frame(Column1 = columns[1], Column2 = columns[2]))
  }
}

# Close the text connection
close(con)

# Print the resulting data frame
print(df)

colnames(df) <- c("motif","tf")

# Function to capitalize the first letter of a string while lower-casing every other letter
# Used to match formatting between mm10 genome and Pando's curated motif list
# Source: https://github.com/quadbio/Pando/issues/30#issue-1564838244 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
  x
}
df$tf <- firstup(df$tf)

# Next get motifs from JASPAR2020
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# For loop to run the workflow through each sample - samples created with 'prepare_pando_objects.RDS'
samples <- c('f.con', 'f.con.c0', 'f.ko', 'f.ko.c0')

# Start the Pando workflow:
for (sam in samples) {
  seurat_object <- readRDS(paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/seurat.objs.RNA.ATAC.', sam, '.RDS'))
  
  # Filter out sequences not in the reference
  RNA <- seurat_object[['RNA']]
  DefaultAssay(seurat_object) <- 'ATAC'
  chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
            "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
  seurat_object <- seurat_object[as.character(seurat_object@assays$ATAC@ranges@seqnames) %in% chrs,]
  seurat_object[['RNA']] <- RNA
  
  # Filter out Rik Genes (these genes start with numbers)
  ATAC <- seurat_object[['ATAC']]
  DefaultAssay(seurat_object) <- 'RNA'
  row_names_with_rik <- rownames(seurat_object@assays$RNA@counts)[grepl("Rik", rownames(seurat_object@assays$RNA@counts))]
  seurat_object <- seurat_object[-(which(rownames(seurat_object@assays$RNA@counts) %in% row_names_with_rik)),]
  seurat_object[['ATAC']] <- ATAC
  
  # # Normalization
  # DefaultAssay(seurat_object) <- 'RNA'
  # seurat_object <- NormalizeData(seurat_object) 
  # seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  # all.genes <- rownames(seurat_object) 
  # seurat_object <- ScaleData(seurat_object, features = all.genes)
  # 
  # DefaultAssay(seurat_object) <- 'ATAC'
  # seurat_object <- RunTFIDF(seurat_object)
  # seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q0')
  # seurat_object <- RunSVD(seurat_object)
  
  # Pando
  seurat_object <- initiate_grn(
    seurat_object,
    rna_assay = 'RNA',
    peak_assay = 'ATAC',
    regions = StringToGRanges(rownames(seurat_object@assays$ATAC))
  )
  regions <- NetworkRegions(seurat_object)
  regions@ranges
  
  seurat_object <- find_motifs(
    seurat_object, 
    genome = BSgenome.Mmusculus.UCSC.mm10,
    pfm = pfm,
    motif_tfs = df
  )
  
  registerDoParallel(8)
  seurat_object <- infer_grn(
    seurat_object,
    peak_to_gene_method = 'GREAT',
    genes = rownames(seurat_object@assays$RNA),
    parallel = T
  )
  
  write.csv(coef(seurat_object), paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/', sam, '_infer_grn.csv'))
  saveRDS(seurat_object, paste0('/storage/chentemp/u250758/mef2c_collab/data/rds/grn/', sam, '.grn.RDS'))
}
