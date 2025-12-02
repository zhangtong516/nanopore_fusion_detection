process MERGE_FASTQ {
    tag { sample }
    storeDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path inpath
    output:
    tuple val(sample), path "*.fastq.gz"
    container 'ubuntu:22.04'
    script:
    """
    set -euo pipefail
    if [ -d "${inpath}" ]; then
      shopt -s nullglob
      files=("${inpath}"/*.fastq "${inpath}"/*.fq "${inpath}"/*.fastq.gz "${inpath}"/*.fq.gz)
      if [ ${#files[@]} -eq 0 ]; then
        echo "No FASTQ files found in ${inpath}" >&2
        exit 1
      fi
      bn=$(basename "${inpath}")
      out="${sample}.${bn}.fastq.gz"
      for f in "${files[@]}"; do
        if [[ "$f" =~ \.gz$ ]]; then
          zcat "$f"
        else
          cat "$f"
        fi
      done | gzip -c > "$out"
    else
      bn=$(basename "${inpath}")
      out="${sample}.${bn}.fastq.gz"
      if [[ "${inpath}" =~ \.gz$ ]]; then
        cp "${inpath}" "$out"
      else
        gzip -c "${inpath}" > "$out"
      fi
    fi
    """
}
