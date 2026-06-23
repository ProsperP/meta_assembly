/*
 * Binning with COMEBin, SemiBin2 and MetaBAT2
 */

def getPowerSet(list) {
    if (list.isEmpty()) return [[]]
    def first = list.head()
    def rest = list.tail()
    def subsets = getPowerSet(rest)
    return subsets + subsets.collect { [first] + it }
}


include { checkm2 } from '../modules/local/checkm/checkm2.nf'
include { summarize_checkm } from '../modules/local/checkm/summarize_checkm.nf'
include { bin_refiner } from '../modules/local/bin_refinement/bin_refiner.nf'


workflow BIN_REFINEMENT {
    take:
    bin_folders
    params

    main:

    //bin_folders.view()
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

    //emit:
    //bins            = ch_bins.ifEmpty([])
    //versions        = ch_versions
}
