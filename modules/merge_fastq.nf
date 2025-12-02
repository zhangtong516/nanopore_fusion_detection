process MERGE_FASTQ {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path inpath
    output:
    tuple val(sample), path "*.fastq.gz"
    container 'ubuntu:22.04'
    script:
    """
    bash ${projectDir}/scripts/merge_fastq.sh "${inpath}" "${sample}"
    """
}
