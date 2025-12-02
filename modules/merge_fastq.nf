process MERGE_FASTQ {
    tag "Merge fastq files"
    storeDir "${params.outdir}/${sample}"
    input:
    tuple val(sample), path (inpath)
    output:
    tuple val(sample), path ("*.fastq.gz")
    script:
    """
    bash ${projectDir}/scripts/merge_fastq.sh "${inpath}" "${sample}"
    """
}
