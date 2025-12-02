process MERGE_FASTQ_FILES {
    tag "Merge all fastq files"
    storeDir "${params.outdir}/${sample}"
    input:
    tuple val(sample), path fastqs
    output:
    tuple val(sample), path "${sample}_merged.fastq.gz"
    container 'ubuntu:22.04'
    script:
    """
    bash ${projectDir}/scripts/merge_fastq_files.sh "${sample}" ${fastqs}
    """
}
