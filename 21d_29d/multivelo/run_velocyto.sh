#!/usr/bin/bash

INPUT=$1

while read -r line
do

SAMPLEID=$line
OUTPUT="/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/"$SAMPLEID
BCFILE="/storage/chentemp/suyangb/m2c_multiome/data/"$SAMPLEID"/filtered_feature_bc_matrix/barcodes.tsv.gz"
BAMFILE="/storage/chentemp/suyangb/m2c_multiome/data/"$SAMPLEID"/gex_possorted_bam.bam"
#BAMFILE="/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/"$SAMPLEID"/cellsorted_possorted_genome_bam.bam"
GENOME="/storage/chentemp/u250758/cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf"

slurmtaco.sh -p short -m 12G -- velocyto run -o $OUTPUT -@ 12 --samtools-memory 10240 --bcfile $BCFILE -e $SAMPLEID $BAMFILE $GENOME

done < $INPUT
