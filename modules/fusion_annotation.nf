process FUSION_ANNOTATION {
    tag { sample }
    publishDir "${params.outdir}/${sample}/", mode: 'move', overwrite: true
    input:
    tuple val(sample), path (fastq)
    output:
    tuple val(sample), path ("${sample}_fusion_sequence.fasta"), path("${sample}_fusion_annotation.gtf")


    script:
    """
    Rscript ${projectDir}/bin/prepare_for_bambu.R ${sample}_jaffa/jaffa_results.csv 
    
    """
}
