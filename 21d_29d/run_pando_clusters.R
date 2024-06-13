library(tidyverse)
library(Pando)
library(Seurat)
library(Signac)
library(doParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)

# Run the pando pipeline on specific clusters from each cluster

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

# For loop to run the workflow through each sample
#samples <- c('21d_m_con', '21d_m_ko', '21d_f_con', '21d_f_ko', '29d_m_con', '29d_m_ko')
samples <- '21d_m_ko'
#clusters <- c('GC0', 'P_Kcnmb2_Calb2', 'P_Ndnf', 'GC_transit_Kcnt2', 'P_Hgf_Stac', 'P_Vwc2l_Reln')
clusters <- c('P_Hgf_Stac_Th', 'P_Vwc2l_Reln_Calb1')

startsWithNumber <- function(string) {
  # Use regular expression to check if the string starts with a number
  grepl("^\\d", string)
}

# Start the Pando workflow:
for (sam in samples) {
  for (clust in clusters) { 
    seurat_object <- readRDS(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam, '.RDS'))
    
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
    seurat_object <- seurat_object[!sapply(rownames(seurat_object@assays$RNA@counts),startsWithNumber),]
    seurat_object[['ATAC']] <- ATAC
    
    seurat_object <- subset(seurat_object, subset = annot.res0.05 == clust)
    
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
      peak_to_gene_method = 'Signac',
      genes = rownames(seurat_object@assays$RNA),
      parallel = T
    )
   
    seurat_object <- find_modules(seurat_object) 
    write.csv(NetworkModules(seurat_object)@meta, paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam, '_', clust, '_grn.csv'))
    saveRDS(seurat_object, paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/', sam, '_', clust, '.grn.RDS'))
    
    #seurat_object <- get_network_graph(seurat_object, features = unique(subset(NetworkModules(seurat_object)@meta, n_genes > 30)$tf)) # Ended up passing in tfs with >30 gene targets
    #png(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/grn_by_cluster/', sam, '_', clust, '.grn.png'))
    #plot_network_graph(seurat_object)
    #dev.off()
  }

    # Start with clusters (GC0, P_Kcnmb2_Calb2, P_Ndnf, P_Hgf_Stac) and TFs (Mef2c Mef2a Bach2)
    seurat_object <- readRDS(paste0(sam, '_', clust, '.grn.RDS'))
    seurat_object <- get_network_graph(seurat_object, graph_name='full_graph', umap_method='none', features = unique(subset(NetworkModules(seurat_object)@meta, n_genes > 10)$tf))
    seurat_object <- get_tf_network(seurat_object, graph='full_graph', tf=tfs[1])
    seurat_object <- get_tf_network(seurat_object, graph='full_graph', tf=tfs[2])
    seurat_object <- get_tf_network(seurat_object, graph='full_graph', tf=tfs[3])
    png(paste0('/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/grn_by_cluster/', unique(seurat_object$orig.ident), '_', unique(seurat_object$annot.res0.05), '_', tfs[1], '.grn.png'))
    plot_tf_network(seurat_object, tf=tfs[1])
    dev.off()
   # Repeat plotting for all three tfs
  }
}