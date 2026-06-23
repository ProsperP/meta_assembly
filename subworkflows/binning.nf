/*
 * Binning with COMEBin, SemiBin2 and MetaBAT2
 */

include { bwa_mem } from '../modules/local/binning/bwa_assembly_align.nf'
include { comebin } from '../modules/local/binning/comebin.nf'
include { semibin2 } from '../modules/local/binning/semibin2.nf'
include { metabat2 } from '../modules/local/binning/metabat2.nf'


workflow BINNING {
    take:
    ch_clean_fqs  // [val{sample}, path(clean_fqs)]
    ch_assembly   // [val{sample}, path(assembly)]
    params

    main:

    ch_versions = channel.empty()
    ch_bins = channel.empty()
    ch_bin_folders = channel.empty()

    ch_bwa_mem_in = channel.empty()
        .mix(ch_clean_fqs).join(ch_assembly)

    bwa_mem(ch_bwa_mem_in)

    ch_binner_in = channel.empty()
        .mix(ch_assembly).join(bwa_mem.out.bam)

    comebin(ch_binner_in)
    ch_versions = ch_versions.mix(comebin.out.versions)

    semibin2(ch_binner_in, params.semibin2)
    ch_versions = ch_versions.mix(semibin2.out.versions)

    metabat2(ch_binner_in, params.metabat2)
    ch_versions = ch_versions.mix(metabat2.out.versions)

    ch_bins = ch_bins.mix(comebin.out.bins)
    ch_bins = ch_bins.join(semibin2.out.bins)
    ch_bins = ch_bins.join(metabat2.out.bins)
    ch_bin_folders = ch_bin_folders
        .mix(comebin.out.bin_folder)
        .mix(semibin2.out.bin_folder)
        .mix(metabat2.out.bin_folder)
        .groupTuple()


    emit:
    assembly_bwa_align = ch_binner_in
    bin_sets          = ch_bins.ifEmpty([])
    bin_folders       = ch_bin_folders
    //unbinned        = ch_unbinned
    //metabat2_depths = ch_combined_depths
    versions        = ch_versions
}
