#!/bin/bash

# Example script to run the nanopore RNA modification pipeline
BASE_DIR="/home/users/astar/gis/zhangt/scratch/nanopore_rna_mod/nanopore_fusion_detection/"
# Set input parameters
INPUT_SAMPLESHEET=$1
OUTPUT_DIR="./jaffa_results_v2.5"

# Create output directory
[[ -d $OUTPUT_DIR ]] || mkdir -p $OUTPUT_DIR

# Run the Nextflow pipeline
nextflow run ${BASE_DIR}/main.nf -c ${BASE_DIR}/nextflow.config \
        --samplesheet $INPUT_SAMPLESHEET --outdir $OUTPUT_DIR \
        -resume

