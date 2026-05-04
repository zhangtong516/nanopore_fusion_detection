process DORADO_ALIGNER {
    tag { sample }
    publishDir "${params.outdir}/${sample}/dorado_align", mode: 'copy'

    input:
    tuple val(sample), path(fastq_file) 
    
    output:
    tuple val(sample), path("${sample}_aligned.bam"), emit: bam

    script:
    """
    dorado aligner \
        ${params.ref_genome} \
        ${fastq_file} \
        --threads ${task.cpus}\
        --mm2-opts '-x splice -k 14' > ${sample}_aligned.bam
    """
}
