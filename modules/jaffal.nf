process RUN_JAFFAL {
    tag { sample }
    publishDir "${params.outdir}/${sample}/jaffal", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path fastq
    output:
    tuple val(sample), dir "${sample}_jaffa"
    container "${params.jaffa_sif ?: 'davidsongroup/jaffa:latest'}"
    script:
    """
    mkdir -p ${sample}_jaffa
    if [ -d "${fastq}" ]; then inputs="${fastq}/*.fastq*"; else inputs="${fastq}"; fi
    bpipe run -d ${sample}_jaffa /JAFFA/JAFFAL.groovy ${inputs}
    """
}
