Pipeline for running Multivelo on Cellranger Multiome output:
1. run_velocyto.sh - First need to run velocyto to get spliced/unspliced rna .loom files
   a. If an error is encountered that suggests sorting the bam files manually using samtools, run the following scripts in order:
     samtools_sort.sh -> rename.sh -> (mv_bam.sh if necessary, moves the sorted renamed bam files to the respective Cellranger outs directory)
   b. Then try run_velocyto.sh again
2. merge_loom.py - merges the created .loom files
3. convert_loom_to_h5ad.py - converts the merged .loom files to an h5ad object (more on h5ad objects on the Anndata website: https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
4. merge_rnaseq.py/merge_atacseq.py - Can be run simultaneously, creates a merged rnaseq and atacseq h5ad object
   a. Run a filtering script like create_gc0_rna_object.py if the pipeline should just be run on a subset of cells (this script subsets just the GC0 cluster)
   b. Requires metadata from another analysis
5. run_multivelo.py - Part 1 of Multivelo
6. seurat_wnn.R - After running part 1, use seurat to find weighted nearest neighbors
7. run_multivelo_pt2.py - Part 2 of Multivelo, uses output from step 5 and 6
8. (Optional) transfer_labels.py - transfer or edit any metadata to the multivelo_results object from another analysis
