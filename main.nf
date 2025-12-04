nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: ''
params.outdir = params.outdir ?: 'results'
params.dorado_model = params.dorado_model ?: 'sup'
params.dorado_device = params.dorado_device ?: 'gpu'
params.jaffa_ref_dir = params.jaffa_ref_dir ?: ''
params.jaffa_sif = params.jaffa_sif ?: ''

include { BASECALL_DORADO } from './modules/dorado'
include { MERGE_FASTQ } from './modules/merge_fastq'
include { MERGE_FASTQ_FILES } from './modules/merge_fastq_files'
include { RUN_JAFFAL } from './modules/jaffal'
include { MAKE_REPORT } from './modules/report'

Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map{ row -> tuple(row.samplename as String, file(row.input_dir)) }
    .set{ samples }

workflow {
    // group multiple rows per sample
    grouped = samples.groupTuple()

    // FASTQ directories per sample
    fastq_dirs = grouped.map{ sample, dirs ->
        def list = dirs.collect{ new java.io.File(it.toString(), 'fastq_pass') }
                        .findAll{ it.exists() && it.isDirectory() && it.listFiles()?.any{ f -> f.name ==~ /(?i).*\.fastq(\.gz)?$/ } }
        tuple(sample, list.collect{ file(it.path) })
    }.filter{ sample, list -> list && list.size() > 0 }

    // POD5 directories per sample, but skip entirely if any fastq_pass exists
    pod5_dirs = grouped.map{ sample, dirs ->
        def hasFastq = dirs.collect{ new java.io.File(it.toString(), 'fastq_pass') }
                           .any{ it.exists() && it.isDirectory() && it.listFiles()?.any{ f -> f.name ==~ /(?i).*\.fastq(\.gz)?$/ } }
        def list = hasFastq ? [] : dirs.collect{ new java.io.File(it.toString(), 'pod5') }
                                  .findAll{ it.exists() && it.isDirectory() }
        tuple(sample, list.collect{ file(it.path) })
    }.filter{ sample, list -> list && list.size() > 0 }

    // Merge FASTQs within each fastq_pass directory -> partial files per sample
    merged_fastq_parts = fastq_dirs.flatMap{ sample, list -> list.collect{ d -> tuple(sample, d) } } | MERGE_FASTQ

    // Basecall POD5 directories -> partial files per sample
    basecalled_parts = pod5_dirs.flatMap{ sample, list -> list.collect{ d -> tuple(sample, d) } } | BASECALL_DORADO

    // Merge partial files per sample into a single FASTQ
    partials = merged_fastq_parts.mix(basecalled_parts)
    final_fastq = partials.groupTuple().map{ sample, items -> tuple(sample, items.collect{ it }) } | MERGE_FASTQ_FILES

    // Run JAFFAL and reporting
    final_fastq | RUN_JAFFAL 
    //| MAKE_REPORT
}
