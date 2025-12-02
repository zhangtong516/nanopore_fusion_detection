process MERGE_FASTQ {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path inpath
    output:
    tuple val(sample), path "${sample}.fastq.gz"
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
      for f in "${files[@]}"; do
        if [[ "$f" =~ \.gz$ ]]; then
          zcat "$f"
        else
          cat "$f"
        fi
      done | gzip -c > ${sample}.fastq.gz
    else
      if [[ "${inpath}" =~ \.gz$ ]]; then
        cp "${inpath}" ${sample}.fastq.gz
      else
        gzip -c "${inpath}" > ${sample}.fastq.gz
      fi
    fi
    """
}

