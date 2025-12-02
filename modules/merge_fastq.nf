process MERGE_FASTQ {
    tag { sample }
    storeDir "${params.outdir}/${sample}"
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
