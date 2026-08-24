/*
 * Binning with COMEBin, SemiBin2, MetaBAT2, MetaDecoder and ConCoCT.
 */

include { bwa_mem } from '../modules/local/binning/bwa_assembly_align'
include { COMEBIN } from '../modules/local/binning/comebin'
include { SEMIBIN2 } from '../modules/local/binning/semibin2'
include { METABAT2 } from '../modules/local/binning/metabat2'
include { METADECODER } from '../modules/local/binning/metadecoder'
include { CONCOCT } from '../modules/local/binning/concoct'
include { CLEAN_BAM } from '../modules/local/binning/clean_bam'
include { collectFlagFiles } from '../modules/local/collect_flag_files'


workflow BINNING {
    take:
    ch_clean_fqs  // [val{sample}, path(clean_fqs)]
    ch_assembly   // [val{sample}, path(assembly)]
    params

    main:

    ch_versions = channel.empty()
    ch_bins = channel.empty()
    ch_bin_folders = channel.empty()

    ch_bwa_mem_in = ch_clean_fqs.join(ch_assembly)
    bwa_mem(ch_bwa_mem_in)

    def binnerCount = 0
    ch_binner_in = ch_assembly.join(bwa_mem.out.bam)

    if ( !params.comebin.skip ) {
        COMEBIN(ch_binner_in)
        ch_versions = ch_versions.mix(COMEBIN.out.versions)
        ch_bins = ch_bins.mix(COMEBIN.out.bins)
        ch_bin_folders = ch_bin_folders.mix(COMEBIN.out.bin_folder)
        binnerCount += 1
    }

    if ( !params.semibin2.skip ) {
        SEMIBIN2(ch_binner_in, params.semibin2)
        ch_versions = ch_versions.mix(SEMIBIN2.out.versions)
        ch_bins = ch_bins.mix(SEMIBIN2.out.bins)
        ch_bin_folders = ch_bin_folders.mix(SEMIBIN2.out.bin_folder)
        binnerCount += 1
    }

    if ( !params.metabat2.skip ) {
        METABAT2(ch_binner_in, params.metabat2)
        ch_versions = ch_versions.mix(METABAT2.out.versions)
        ch_bins = ch_bins.mix(METABAT2.out.bins)
        ch_bin_folders = ch_bin_folders.mix(METABAT2.out.bin_folder)
        binnerCount += 1
    }

    if ( !params.metadecoder.skip ) {
        METADECODER(ch_binner_in, params.metadecoder)
        ch_versions = ch_versions.mix(METADECODER.out.versions)
        ch_bins = ch_bins.mix(METADECODER.out.bins)
        ch_bin_folders = ch_bin_folders.mix(METADECODER.out.bin_folder)
        binnerCount += 1
    }

    if ( !params.concoct.skip ) {
        CONCOCT(ch_binner_in, params.concoct)
        ch_versions = ch_versions.mix(CONCOCT.out.versions)
        ch_bins = ch_bins.mix(CONCOCT.out.bins)
        ch_bin_folders = ch_bin_folders.mix(CONCOCT.out.bin_folder)
        binnerCount += 1
    }

    log.info "Number of enabled binners: ${binnerCount}"

    ch_bins = ch_bins.groupTuple(size: binnerCount)
    ch_bin_folders = ch_bin_folders.groupTuple(size: binnerCount)

    CLEAN_BAM(ch_bin_folders.join(bwa_mem.out.bam).join(bwa_mem.out.flag_file))

    emit:
    //assembly_bwa_align = ch_binner_in
    bin_sets          = ch_bins.ifEmpty([])
    bin_folders       = ch_bin_folders
    //unbinned        = ch_unbinned
    //metabat2_depths = ch_combined_depths
    versions        = ch_versions
}