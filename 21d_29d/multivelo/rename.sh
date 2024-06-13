#!/usr/bin/bash

INPUT=$1

while read -r line
do

SAMPLEID=$line
OUTPUT="/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/"$SAMPLEID"/cellsorted_gex_possorted_bam.bam"
BAMFILE="/storage/chentemp/u250758/mef2c_collab/data/21d_29d_all/loom/"$SAMPLEID"/cellsorted_possorted_genome_bam.bam"

slurmtaco.sh -p short -- mv $BAMFILE $OUTPUT

done < $INPUT
