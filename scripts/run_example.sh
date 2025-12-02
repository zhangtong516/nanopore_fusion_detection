#!/usr/bin/env bash
set -euo pipefail
SAMPLESHEET=${SAMPLESHEET:-samples.csv}
OUTDIR=${OUTDIR:-results}
DORADO_MODEL=${DORADO_MODEL:-/path/to/dorado/model}
DORADO_DEVICE=${DORADO_DEVICE:-gpu}
JAFFA_REF_DIR=${JAFFA_REF_DIR:-/path/to/jaffa_refs}
JAFFA_SIF=${JAFFA_SIF:-/path/to/jaffa.sif}
nextflow run main.nf -profile apptainer --samplesheet "$SAMPLESHEET" --outdir "$OUTDIR" --dorado_model "$DORADO_MODEL" --dorado_device "$DORADO_DEVICE" --jaffa_ref_dir "$JAFFA_REF_DIR" --jaffa_sif "$JAFFA_SIF"

