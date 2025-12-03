process MAKE_REPORT {
    tag "Summary of fusion detection"
    publishDir "${params.outdir}/${sample}", mode="copy" 
    input:
    tuple val(sample), path (jaffal_dir)
    output:
    path "${sample}_report.html"
    script:
    """
    python ${projectDir}/bin/make_report.py --sample ${sample} --jaffal_dir ${jaffal_dir} --out ${sample}_report.html
    """
}
