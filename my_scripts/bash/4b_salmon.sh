#!/bin/bash

INDEXES=( "Ref" "Gamma" "Proton" )
INPUT_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022/akinshin_sd"
INDEX_DIR="/home/akinshin_sd/hvulgare/reference/transcriptome/hvulgare_transcriptome_index/"
OUTPUT_DIR="/home/akinshin_sd/hvulgare/salmon"

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  R1_INPUT="$INPUT_DIR/${FULL_NAME}_trimR1.fastq"
  R2_INPUT="$INPUT_DIR/${FULL_NAME}_trimR2.fastq"
  OUTPUT="$OUTPUT_DIR/${FULL_NAME}.sf"

  salmon quant -i $INDEX_DIR -l A -1 ${R1_INPUT} -2 ${R2_INPUT} -o ${OUTPUT}
  mv "${OUTPUT}/quant.sf" "${OUTPUT}/${FULL_NAME}_quant.sf"
  done
done
