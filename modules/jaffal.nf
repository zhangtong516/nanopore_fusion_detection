process RUN_JAFFAL {
    tag { sample }
    storeDir "${params.outdir}/${sample}/jaffal", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path fastq
    output:
    tuple val(sample), dir "${sample}_jaffa"
    script:
    """
    set -euo pipefail
    mkdir -p ${sample}_jaffa
    cd ${sample}_jaffa
    img="${params.jaffa_sif ?: 'docker://davidsongroup/jaffa:latest'}"
    bind_ref="${params.jaffa_ref_dir ?: ''}"
    fastq_dir=$(dirname "${fastq}")
    BREF=""; [ -n "$bind_ref" ] && BREF="-B ${bind_ref}:/ref"
    apptainer run $BREF -B "$fastq_dir":"$fastq_dir" -B "$(pwd)":"$(pwd)" "$img" /JAFFA/JAFFAL.groovy "${fastq}"
    """
}
