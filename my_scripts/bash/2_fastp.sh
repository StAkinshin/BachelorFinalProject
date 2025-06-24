#!/bin/bash

DATA_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022"

INDEXES=( "Ref" "Gamma" "Proton" )
OUTPUT_DIR="/data/illumina_seq/transcriptome/hordeum_vulgare_2022/akinshin_sd"
REPORT_DIR="/home/akinshin_sd/hvulgare/fastp"
AVERAGE_QUAL=30
THREAD=16

for index in "${INDEXES[@]}"; do

FULL_NAME="HV_Rep3_${index}"
R1_INPUT="$DATA_DIR/${FULL_NAME}_R1.fastq"
R2_INPUT="$DATA_DIR/${FULL_NAME}_R2.fastq"
R1_OUTPUT="$REPORT_DIR/${FULL_NAME}_trimR1.fastq"
R2_OUTPUT="$REPORT_DIR/${FULL_NAME}_trimR2.fastq"
REPORT_NAME="$REPORT_DIR/${FULL_NAME}_fastq.html"

fastp --thread $THREAD --in1 ${R1_INPUT} --out1 ${R1_OUTPUT} --in2 ${R2_INPUT} --out2 ${R2_OUTPUT} --average_qual $AVERAGE_QUAL --html ${REPORT_NAME} 
done
 
