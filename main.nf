#!/usr/bin/env nextflow


// include pre-defined modules
include { kneaddata } from './modules/local/kneaddata.nf'
include { QUAST } from './modules/local/quast'

include { ASSEMBLY } from './subworkflows/assembly'
include { BINNING } from './subworkflows/binning.nf'
include { BIN_REFINEMENT } from './subworkflows/bin_refinement.nf'


// define user-supplied parameters
params {
    samples: Path
    kneaddata: Map
    spades: Map
    quast: Map
    binning: Map
}

workflow {

    main:
    ch_clean_reads = channel.empty()
    ch_kneaddata_log = channel.empty()

    input_fastqs_ch = channel.fromPath(params.samples)
        .splitCsv(sep: '\t')
        .groupBy()
        .map { key, vals -> tuple(key, vals.toSorted()) }


    if ( params.kneaddata.skip ) {
        ch_clean_reads = input_fastqs_ch
    } else {
        kneaddata(input_fastqs_ch, params.kneaddata)
        ch_clean_reads = ch_clean_reads.mix(kneaddata.out.clean_fqs)
        ch_kneaddata_log = ch_kneaddata_log.mix(kneaddata.out.log)
    }

    ASSEMBLY(ch_clean_reads, params.assembly)
    ch_final_assembly = ASSEMBLY.out.final_assembly
    ch_assembly_out = ASSEMBLY.out.assembly_out
    ch_assembly_log = ASSEMBLY.out.assembly_log

    QUAST(ch_final_assembly)
    ch_quast_report = QUAST.out.assembly_report

    BINNING(ch_clean_reads, ch_final_assembly, params.binning)
    BIN_REFINEMENT(
        ch_final_assembly,
        BINNING.out.bin_folders,
        params.binning
    )

    publish:
    clean_fqs = ch_clean_reads.ifEmpty([])
    kneaddata_log = ch_kneaddata_log.ifEmpty([])
    assembly_out = ch_assembly_out
    assembly_log = ch_assembly_log
    final_assembly = ch_final_assembly
    quast_report = ch_quast_report
    assembly_bwa_align = BINNING.out.assembly_bwa_align
    bin_sets = BINNING.out.bin_sets
    //versions = BINNING.out.versionsa

}

output {
    clean_fqs { path "${params.kneaddata.outdir}/" }
    kneaddata_log { path "${params.kneaddata.outdir}/" }
    assembly_out { path "${params.assembly.outdir}" }
    assembly_log { path "${params.assembly.outdir}" }
    final_assembly { path "${params.assembly.outdir}" }
    quast_report { path "${params.assembly.outdir}" }
    assembly_bwa_align { path "${params.binning.outdir}" }
    bin_sets { path "${params.binning.outdir}" }
    //versions { path "${params.binning.outdir}" }
}
