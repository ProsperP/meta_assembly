/*
 * Bin refinement using metaWRAP or MAGScoT
 */

def getPowerSet(list) {
    if (list.isEmpty()) return [[]]
    def first = list.head()
    def rest = list.tail()
    def subsets = getPowerSet(rest)
    return subsets + subsets.collect { it -> [first] + it }
}


include { identify_single_copy_marker } from '../modules/local/bin_refinement/identify_single_copy_marker.nf'
include {
    format_contigs2bin
    magscot
} from '../modules/local/bin_refinement/magscot.nf'
include { checkm2 } from '../modules/local/checkm/checkm2.nf'
include { summarize_checkm } from '../modules/local/checkm/summarize_checkm.nf'
include { bin_refiner } from '../modules/local/bin_refinement/bin_refiner.nf'


workflow BIN_REFINEMENT {
    take:
    assembly
    bin_folders
    params

    main:

    def binnerCount = 0
    if ( !params.comebin.skip ) binnerCount += 1
    if ( !params.semibin2.skip ) binnerCount += 1
    if ( !params.metabat2.skip ) binnerCount += 1
    if ( !params.metadecoder.skip ) binnerCount += 1
    if ( !params.concoct.skip ) binnerCount += 1

    //bin_folders.view()
    ch_refined_bins = channel.empty()
    ch_refine_stats = channel.empty()
    ch_versions = channel.empty()

    if( params.refinement.method == 'metaWRAP' ) {
        ch_in_bin_refiner = bin_folders.map { sample, binners, bins_dirs ->
                def bins_map = [binners, bins_dirs].transpose().collectEntries()
                return [sample, bins_map]
            }
            .flatMap { sample, bins_map ->
                def binner_list = bins_map.keySet().toList()
                def all_subsets = getPowerSet(binner_list)
                def combinations = []
            
                // loop through all subsets, keep only set size larger than 1
                // since bins from a single binner does not need to be merged
                all_subsets.each { subset ->
                    if (subset.size() >= 2) {
                        def combo_name = subset.join('_')
                        // Extract corresponding bin path from the binner
                        def combo_dirs = subset.collect { binner -> bins_map[binner] }

                        combinations << [sample, combo_name, combo_dirs]
                    }
                }
                return combinations // each element in the collection is emitted separately
            }
        bin_refiner(ch_in_bin_refiner)

        ch_in_bin_checkm2 = bin_folders.flatMap { sample, binners, bins_dirs ->
                [binners, bins_dirs].transpose()
                    .collect { pair -> [sample] + pair }
            }
            .mix(bin_refiner.out.bins_comb)

        checkm2(ch_in_bin_checkm2, params.checkm2)
        summarize_checkm(checkm2.out.checkm2_tsv)
        //ch_bins_comb_checkm2 = checkm2(bin_refiner.out.bins_comb)
    } else if( params.refinement.method == 'MAGScoT' ) {
        identify_single_copy_marker(assembly, params.refinement.magscot)
        ch_versions = ch_versions.mix(identify_single_copy_marker.out.versions)

        ch_in_format_contig2bin = bin_folders.transpose()
        format_contigs2bin(ch_in_format_contig2bin)

        ch_in_magscot = assembly.join(
            format_contigs2bin.out
                .groupTuple(size: binnerCount)
                .map { sample, _binners, bins_dirs -> [sample, bins_dirs] }
            )
            .join(identify_single_copy_marker.out.hmm)

        magscot(ch_in_magscot, params.refinement.magscot)
        ch_refined_bins = ch_refined_bins.mix(magscot.out.refined_bins)
        ch_refine_stats = ch_refine_stats.mix(magscot.out.magscot_out)
    }

    emit:
    refined_bins = ch_refined_bins.ifEmpty([])
    refine_stats = ch_refine_stats.ifEmpty([])
    versions     = ch_versions
}