process RUN_JAFFAL {
    tag { sample }
    publishDir "${params.outdir}/${sample}/jaffal", mode: 'move', overwrite: true
    input:
    tuple val(sample), path (fastq)
    output:
    tuple val(sample), path ("${sample}_jaffa")


    script:
    """
    set -euo pipefail
    mkdir -p ${sample}_jaffa
    cd ${sample}_jaffa

    apptainer run -B ${params.jaffa_ref_dir}:/ref/ ${params.jaffa_sif} /JAFFA/JAFFAL.groovy ../${fastq}

    """
}

