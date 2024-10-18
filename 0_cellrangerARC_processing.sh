#!/bin/bash
# example usage ./0_cellrangerARC_processing.sh 1dpa1
# Check if a sample ID was provided as an input argument
if [ -z "$1" ]; then
  echo "Please provide a sample ID as an argument."
  exit 1
fi

# Assign the sample ID to a variable
SAMPLE_ID="$1"

# Setting PATH to include CellRanger-ARC directory
export PATH="/bar/ichen/cellranger-arc-2.0.0/:$PATH"

# Running cellranger-arc with appropriate arguments
nohup cellranger-arc count --id=$SAMPLE_ID \
  --reference=refdata-cellranger-arc-2.0.0-danRer11_onlyChr_wEGFP/ \
  --libraries=${SAMPLE_ID}_10xMultiome_lib.csv \
  --localcores=23 \
  --localmem=90 \
  > ${SAMPLE_ID}_cellranger-arc-2.0.0.log 2>&1 &
