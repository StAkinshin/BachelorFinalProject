#!/bin/bash

INDEXES=( "Ref" "Gamma" "Proton" )
INPUT_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022/akinshin_sd"
INDEX="/home/akinshin_sd/hvulgare/reference/transcriptome/hvulgare_transcriptome.idx"
OUTPUT_DIR="/home/akinshin_sd/hvulgare/kallisto"
THREAD=24

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  R1_INPUT="$INPUT_DIR/${FULL_NAME}_trimR1.fastq"
  R2_INPUT="$INPUT_DIR/${FULL_NAME}_trimR2.fastq"
  OUTPUT="$OUTPUT_DIR/${FULL_NAME}"

  kallisto quant -t $THREAD -i $INDEX -o ${OUTPUT} -b 100 ${R1_INPUT} ${R2_INPUT}
  done
done
