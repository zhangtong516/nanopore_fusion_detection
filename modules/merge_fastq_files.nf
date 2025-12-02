process MERGE_FASTQ_FILES {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path (fastqs)
    output:
    tuple val(sample), path("${sample}_merged.fastq.gz")
    script:
    """
    set -euo pipefail
    if [ -z "${fastqs}" ]; then echo "No FASTQ inputs" >&2; exit 1; fi
    for f in ${fastqs}; do
      if [[ "\$f" =~ \\.gz\$ ]]; then zcat "\$f"; else cat "\$f"; fi
    done | gzip -c > ${sample}_merged.fastq.gz
    """
}
