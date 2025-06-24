#!/bin/bash

INDEXES=( "Ref" "Gamma" "Proton" )
DATA_DIR="/home/akinshin_sd/hvulgare/featureCounts"

THREAD=16

for index in "${INDEXES[@]}"; do
  for rep in {1..3}; do
  FULL_NAME="HV_Rep${rep}_${index}"
  INPUT_DATA="$DATA_DIR/${FULL_NAME}_counts.txt"
  OUTPUT_DATA="$DATA_DIR/${FULL_NAME}_cutCounts.txt"
  cat ${INPUT_DATA} | cut -f 1,7 | less > ${OUTPUT_DATA}
  done
done
 
