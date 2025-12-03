#!/usr/bin/env bash
set -euo pipefail


nextflow run main.nf \
     --samplesheet "$SAMPLESHEET" \
     --outdir "$OUTDIR" \
     --profile apptainer \
     -resume 

