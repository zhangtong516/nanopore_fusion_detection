process BASECALL_DORADO {
    tag { sample }
    publishDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path pod5_dir
    output:
    tuple val(sample), path "${sample}.fastq.gz"
    container 'ontresearch/dorado:latest'
    script:
    """
    dorado basecaller --device ${params.dorado_device} ${params.dorado_model} ${pod5_dir} | gzip -c > ${sample}.fastq.gz
    """
}

