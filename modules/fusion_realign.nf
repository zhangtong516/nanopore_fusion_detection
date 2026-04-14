process FUSION_REALIGN {
    tag { sample }
    publishDir "${params.outdir}/${sample}/", mode: 'move', overwrite: true
    input:
    tuple val(sample), path(bam), path (fusion_fasta), path(fusion_bed)
    output:
    tuple val(sample), path ("${sample}_fusion_sequence.fasta"), path("${sample}_fusion_annotation.gtf")


    script:
    """
    samtools view -bhL ${fusion_bed} ${bam} | samtools fastq - > ${sample}_fusion.fastq
    minimap2 -ax splice -G2200k -N 5 --sam-hit-only -t 20 ${fusion_fasta} ${sample}_fusion.fastq | samtools sort -O bam -o ${sample}_fusion.bam -
        samtools index ${sample}_fusion.bam
        
    """
}
