# Nanopore Fusion Detection Pipeline

This Nextflow pipeline detects fusion genes from Nanopore long-read transcriptome data using JAFFAL. It supports two input modes per sample: using existing `fastq_pass` reads directly, or basecalling `pod5` with Dorado. All FASTQ files are merged into a single `sample.fastq.gz` before running JAFFAL.

## Sample Sheet

- CSV headers: `sample,sample_dir`
- `sample_dir` must be a folder that may contain `fastq_pass` and/or `pod5` subfolders.
- Example:

```
sample,sample_dir
S1,/data/samples/S1
S2,/data/samples/S2
```

## Run (Apptainer)

```
nextflow run main.nf -profile apptainer \
  --samplesheet samples.csv \
  --outdir results \
  --dorado_model /path/to/dorado/model \
  --dorado_device gpu \
  --jaffa_ref_dir /path/to/jaffa_refs \
  --jaffa_sif /path/to/jaffa.sif
```

Or use the provided script:

```
bash scripts/run_example.sh
```

## Run (Docker)

```
nextflow run main.nf -profile docker \
  --samplesheet samples.csv \
  --outdir results \
  --dorado_model /path/to/dorado/model \
  --dorado_device gpu \
  --jaffa_ref_dir /path/to/jaffa_refs
```

## Outputs

- `results/<sample>/<sample>.fastq.gz`: merged or basecalled FASTQ.
- `results/<sample>/jaffal/`: JAFFAL outputs.
- `results/<sample>/<sample>_report.html`: HTML summary with links and previews.

## Notes

- Dorado requires a valid model path (`--dorado_model`).
- JAFFAL references must be accessible inside the container at `/ref` (`--jaffa_ref_dir` binds host refs to `/ref`).
- The pipeline auto-detects `fastq_pass`; if not present, it uses `pod5` for basecalling. If neither exists, the sample errors.
