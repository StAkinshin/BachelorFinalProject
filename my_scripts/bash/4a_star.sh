#!/bin/bash

DATA_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022/akinshin_sd"

INDEXES=( "Ref" "Gamma" "Proton" )
OUTPUT_DIR="/home/akinshin_sd/hvulgare/star"
GENOME_INDEX="/home/akinshin_sd/hvulgare/reference/genome"
THREAD=16

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  R1_DATA="$DATA_DIR/${FULL_NAME}_trimR1.fastq"
  R2_DATA="$DATA_DIR/${FULL_NAME}_trimR2.fastq"
  STAR --runThreadN $THREAD --genomeDir $GENOME_INDEX --readFilesIn ${R1_DATA} ${R2_DATA} --outFileNamePrefix $OUTPUT_DIR/${FULL_NAME}_ --outSAMtype BAM SortedByCoordinate 
  done
done
 
