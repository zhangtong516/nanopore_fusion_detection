process DORADO_ALIGNER {
    tag { sample }
    publishDir "${params.outdir}/${sample}/dorado_align", mode: 'copy'

    input:
    tuple val(sample), path(fastq_file) 
    output:
    tuple val(sample), path("${samplen}_aligned.bam"), emit: bam

    script:
    """
    dorado aligner \
        ${reference_genome} \
        ${fastq_file} \
        --threads ${params.jaffa_ref_dir}/hg38.fa \
        --mm2-opts '-x splice -k 14' > ${samplen}_aligned.bam
    """
}
