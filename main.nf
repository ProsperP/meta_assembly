#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include pre-defined modules
include { kneaddata } from './modules/local/kneaddata.nf'
include { FASTQ_SORT } from './modules/local/fastq_sort'
include { QUAST } from './modules/local/quast'

include { ASSEMBLY } from './subworkflows/assembly'
include { BINNING } from './subworkflows/binning.nf'
include { BIN_REFINEMENT } from './subworkflows/bin_refinement.nf'
include { DEREPLICATION } from './subworkflows/dereplication'
include { QUANT_BINS } from './subworkflows/quant_bins'
include { GTDBTK_CLASSIFY_WF } from './modules/local/gtdbtk'
include { IQTREE3 } from './modules/local/iqtree'
include { BAKTA } from './modules/local/bakta'


// define user-supplied parameters
params {
    samples: Path
    kneaddata: Map
    spades: Map
    quast: Map
    binning: Map
    dereplication: Map
    annotation: Map
    drep_batch_size: Integer = 500
    funct_batch_size: Integer = 50
    coverm: Map
}

workflow {

    main:
    ch_clean_reads = channel.empty()
    ch_kneaddata_log = channel.empty()
    ch_assembly_out = channel.empty()
    ch_assembly_log = channel.empty()

    input_fastqs_ch = channel.fromPath(params.samples)
        .splitCsv(sep: '\t')
        .groupBy()
        .map { key, vals -> tuple(key, vals.toSorted()) }


    if ( params.kneaddata.skip ) {
        ch_clean_reads = FASTQ_SORT(input_fastqs_ch)
    } else {
        kneaddata(input_fastqs_ch, params.kneaddata)
        FASTQ_SORT(kneaddata.out.clean_fqs)
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

    DEREPLICATION(
        BIN_REFINEMENT.out.refined_bins,
        params.dereplication,
        params.drep_batch_size
    )

    QUANT_BINS(DEREPLICATION.out.mags, ch_clean_reads, params.coverm)

    GTDBTK_CLASSIFY_WF(DEREPLICATION.out.mags)
    IQTREE3(GTDBTK_CLASSIFY_WF.out.user_msa)
    ch_bakta_in = DEREPLICATION.out.mags
        .flatten()
        .collate( params.funct_batch_size )
    BAKTA(ch_bakta_in, params.annotation.bakta)

    publish:
    clean_fqs = ch_clean_reads.ifEmpty([])
    kneaddata_log = ch_kneaddata_log.ifEmpty([])
    assembly_out = ch_assembly_out
    assembly_log = ch_assembly_log
    final_assembly = ch_final_assembly
    quast_report = ch_quast_report
    bin_sets = BINNING.out.bin_sets
    refined_bins = BIN_REFINEMENT.out.refined_bins
    refine_stats = BIN_REFINEMENT.out.refine_stats
    mags = DEREPLICATION.out.mags.mix(DEREPLICATION.out.id_convert)
    srgs = DEREPLICATION.out.srgs
    mags_profile = QUANT_BINS.out
    taxonomy = GTDBTK_CLASSIFY_WF.out.taxonomy
    tree = IQTREE3.out
    function = BAKTA.out.collect()

}

output {
    clean_fqs { path "${params.kneaddata.outdir}/" }
    kneaddata_log { path "${params.kneaddata.outdir}/" }
    assembly_out { path "${params.assembly.outdir}" }
    assembly_log { path "${params.assembly.outdir}" }
    final_assembly { path "${params.assembly.outdir}" }
    quast_report { path "${params.assembly.outdir}" }
    bin_sets { path "${params.binning.outdir}" }
    refined_bins { path "${params.binning.outdir}" }
    refine_stats { path "${params.binning.outdir}" }
    mags { path "${params.binning.outdir}" }
    srgs { path "${params.binning.outdir}/SRGs" }
    mags_profile { path "${params.coverm.outdir}" }
    taxonomy { path "${params.annotation.outdir}/taxonomy" }
    tree { path "${params.annotation.outdir}" }
    function { path "${params.annotation.outdir}/function" }
}
