#!/usr/bin/env nextflow


// include pre-defined modules
include { kneaddata } from './modules/local/kneaddata.nf'
include {
    spades
    filter_assembly
    quast
} from './modules/local/assembly.nf'
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

    spades(ch_clean_reads, params.spades)
    filter_assembly(spades.out.scaffolds, params.spades)
    quast(filter_assembly.out.final_assembly)
    BINNING(ch_clean_reads, filter_assembly.out.final_assembly, params.binning)
    BIN_REFINEMENT(BINNING.out.bin_folders, params.binning)

    publish:
    clean_fqs = ch_clean_reads.ifEmpty([])
    kneaddata_log = ch_kneaddata_log.ifEmpty([])
    spades_out = channel.empty().mix(
        spades.out.contigs,
        spades.out.scaffolds,
        spades.out.gene_clusters,
        spades.out.transcripts,
        spades.out.gfa
    ).ifEmpty([])
    spades_logs = channel.empty().mix(
        spades.out.warnings,
        spades.out.log
    )
    final_assembly = filter_assembly.out.final_assembly
    quast_report = quast.out.assembly_report
    assembly_bwa_align = BINNING.out.assembly_bwa_align
    bin_sets = BINNING.out.bin_sets
    //versions = BINNING.out.versionsa

}

output {
    clean_fqs { path "${params.kneaddata.outdir}/" }
    kneaddata_log { path "${params.kneaddata.outdir}/" }
    spades_out { path "${params.spades.outdir}" }
    spades_logs { path "${params.spades.outdir}" }
    final_assembly { path "${params.spades.outdir}" }
    quast_report { path "${params.spades.outdir}" }
    assembly_bwa_align { path "${params.binning.outdir}" }
    bin_sets { path "${params.binning.outdir}" }
    //versions { path "${params.binning.outdir}" }
}
