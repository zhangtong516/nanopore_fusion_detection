process MERGE_FASTQ {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path(inpath)
    output:
    tuple val(sample), path("*.fastq.gz")
    script:
    """
    bn=\$(basename ${inpath}) 
    find ./${inpath}/ -name "*fastq.gz" -exec cat {} + > ${sample}_\${bn}.fastq.gz 
    """
}
