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
    img="${params.jaffa_sif}"
    BREF=""; [ -n "${params.jaffa_ref_dir}" ] && BREF="-B ${params.jaffa_ref_dir}:/ref/"
    fqdir=$(dirname "${fastq}")
    apptainer run $BREF -B "$fqdir":"$fqdir" -B "$(pwd)":"$(pwd)" "$img" /JAFFA/JAFFAL.groovy "${fastq}"
    """
}
