#!/bin/bash

INDEXES=( "Ref" "Gamma" "Proton" )
INPUT_DIR="/home/akinshin_sd/hvulgare/star"
GENOME_GTF="/home/akinshin_sd/hvulgare/reference/genome/hvulgare_genome.gtf"
OUTPUT_DIR="/home/akinshin_sd/hvulgare/featureCounts"

THREAD=16

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  INPUT_DATA="$INPUT_DIR/${FULL_NAME}_aligned.bam"
  OUTPUT_DATA="$OUTPUT_DIR/${FULL_NAME}_counts.txt"
  featureCounts -T $THREAD -a $GENOME_GTF -p -o ${OUTPUT_DATA} --largestOverlap -t exon -g gene_id ${INPUT_DATA} 
  done
done
 
