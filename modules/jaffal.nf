process RUN_JAFFAL {
    tag "Run JAFFA for fusion"
    storeDir "${params.outdir}/${sample}/"
    input:
    tuple val(sample), path fastq
    output:
    tuple val(sample), dir "${sample}_jaffa"
    script:
    """
    set -euo pipefail
    mkdir -p ${sample}_jaffa
    cd ${sample}_jaffa
    apptainer run $BREF -B ${params.jaffa_ref_dir}:/ref/  ${params.jaffa_sif} /JAFFA/JAFFAL.groovy ${fastq}
    """
}

// BREF=""; [ -n "$bind_ref" ] && BREF="-B ${bind_ref}:/ref"