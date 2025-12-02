nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: ''
params.outdir = params.outdir ?: 'results'
params.dorado_model = params.dorado_model ?: ''
params.dorado_device = params.dorado_device ?: 'gpu'
params.jaffa_ref_dir = params.jaffa_ref_dir ?: ''
params.jaffa_sif = params.jaffa_sif ?: ''

include { BASECALL_DORADO } from './modules/dorado.nf'
include { MERGE_FASTQ } from './modules/merge_fastq.nf'
include { RUN_JAFFAL } from './modules/jaffal.nf'
include { MAKE_REPORT } from './modules/report.nf'

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map{ row -> tuple(row.sample as String, file(row.sample_dir)) }
    .set{ samples }

workflow {
    samples.into { samples_fastq; samples_pod5 }
    samples_fastq = samples_fastq.filter{ sample, dir ->
        def fq = new java.io.File(dir.toString(), 'fastq_pass')
        fq.exists() && fq.isDirectory() && fq.listFiles()?.any{ f -> f.name ==~ /(?i).*\.fastq(\.gz)?$/ }
    }.map{ sample, dir -> tuple(sample, file("${dir}/fastq_pass")) }

    samples_pod5 = samples_pod5.filter{ sample, dir ->
        def fq = new java.io.File(dir.toString(), 'fastq_pass')
        def p5 = new java.io.File(dir.toString(), 'pod5')
        def hasFastq = fq.exists() && fq.isDirectory() && fq.listFiles()?.any{ f -> f.name ==~ /(?i).*\.fastq(\.gz)?$/ }
        def hasPod5 = p5.exists() && p5.isDirectory()
        if (!hasFastq && !hasPod5)
            throw new IllegalArgumentException("Sample ${sample}: neither fastq_pass nor pod5 folder found in ${dir}")
        return !hasFastq && hasPod5
    }.map{ sample, dir -> tuple(sample, file("${dir}/pod5")) }

    basecalled = samples_pod5 | BASECALL_DORADO
    merged_from_dir = samples_fastq | MERGE_FASTQ
    merged_all = merged_from_dir.mix(basecalled)
    merged_all | RUN_JAFFAL | MAKE_REPORT
}

