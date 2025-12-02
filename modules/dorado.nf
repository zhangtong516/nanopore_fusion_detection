process BASECALL_DORADO {
    tag "Basecall with Dorado"
    storeDir "${params.outdir}/${sample}" 
    input:
    tuple val(sample), path pod5_dir
    output:
    tuple val(sample), path "*.fastq.gz"
    container 'ontresearch/dorado:latest'
    script:
    """
    bn=`basename "${pod5_dir}"`
    dorado basecaller --device ${params.dorado_device} ${params.dorado_model} ${pod5_dir} | gzip -c > ${sample}.${bn}.fastq.gz
    """
}
