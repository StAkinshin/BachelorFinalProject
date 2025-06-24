#!/bin/bash

DATA_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022/akinshin_sd"

INDEXES=( "Ref" "Gamma" "Proton" )
OUTPUT_DIR="/home/akinshin_sd/hvulgare/fastqc/fastp_qc"

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  R1_DATA="$DATA_DIR/${FULL_NAME}_trimR1.fastq"
  R2_DATA="$DATA_DIR/${FULL_NAME}_trimR2.fastq"
  fastqc --threads 16 --outdir $OUTPUT_DIR ${R1_DATA} ${R2_DATA}
  done
done
 
