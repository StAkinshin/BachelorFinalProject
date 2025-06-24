#!/bin/bash

DATA_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022"

INDEXES=( "Ref" "Gamma" "Proton" )
OUTPUT_DIR="/home/akinshin_sd/hvulgare/fastqc/raw_qc"

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  R1_DATA="$DATA_DIR/${FULL_NAME}_R1.fastq"
  R2_DATA="$DATA_DIR/${FULL_NAME}_R2.fastq"
  fastqc --threads 16 --outdir $OUTPUT_DIR ${R1_DATA} ${R2_DATA}
  done
done
 
