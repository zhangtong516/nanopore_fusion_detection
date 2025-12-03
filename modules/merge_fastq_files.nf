process MERGE_FASTQ_FILES {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'symlink', overwrite: true
    input:
    tuple val(sample), path(fastqs)
    output:
    tuple val(sample), path("${sample}_merged.fastq.gz")
    script:
    """
    set -euo pipefail
    files=( ${fastqs} )
    if [ \${#files[@]} -eq 0 ]; then
      echo "No FASTQ inputs" >&2
      exit 1
    elif [ \${#files[@]} -eq 1 ]; then
      f="\${files[0]}"
      cp "\$f" ${sample}_merged.fastq.gz
    else
      for f in "\${files[@]}"; do
        cat "\$f"
      done > ${sample}_merged.fastq.gz
    fi
    """
