process FUSION_ANNOTATION {
    tag { sample }
    publishDir "${params.outdir}/${sample}/fusion_annotation/", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(jaffa_dir), path(bam)

    
    output:
    tuple val(sample), path("${sample}_bambu.html"), path("${sample}_fusionGene_sequences.fa"), path("${sample}_fusion.bam"), path("${sample}_fusion.bam.bai"), emit: annotation
    

    script:
    """
    set -euo pipefail
    
    # Check if GTF and reference genome are provided
    if [[ -z "${params.gtf}" ]] || [[ -z "${params.ref_genome}" ]]; then
        echo "Error: params.gtf and params.ref_genome must be specified for fusion annotation"
        exit 1
    fi

    apptainer run \
        -B "${PWD}:${PWD}" \
        -B "${projectDir}:${projectDir}" \
        -B "${params.gtf}:${params.gtf}" \
        -B "${params.ref_genome}:${params.ref_genome}" \
        ${params.bambu_sif} Rscript \
        "${projectDir}/scripts/render_bambu_report.R" \
        --samplename "${sample}" \
        --bamfile "${bam}" \
        --jaffaresults "${jaffa_dir}/jaffa_results.csv" \
        --gtf "${params.gtf}" \
        --refgenome "${params.ref_genome}" \
        --rmdfile "${projectDir}/bin/bambu.Rmd" 
    """
}


