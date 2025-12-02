process MAKE_REPORT {
    tag "Summary of fusion detection"
    storeDir "${params.outdir}/${sample}", mode: 'copy', overwrite: true
    input:
    tuple val(sample), path jaffal_dir
    output:
    path "${sample}_report.html"
    script:
    """
    python ${projectDir}/bin/make_report.py --sample ${sample} --jaffal_dir ${jaffal_dir} --out ${sample}_report.html
    """
}
