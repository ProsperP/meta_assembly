/*
 Quant bins profile using coverm
 */


include { CATENATE_BINS } from '../modules/local/catenate_bins'
include { BWA_INDEX } from '../modules/local/bwa/bwa_index'
include { BWA_MEM } from '../modules/local/bwa/bwa_mem'
include {
    COVERM
    MERGE_COVERM
} from '../modules/local/coverm'
include { CLEAN_BAM } from '../modules/local/bwa/clean_bam.nf'


workflow QUANT_BINS {
    take:
    bins
    clean_fqs
    params

    main:

    CATENATE_BINS(bins)
    BWA_INDEX(CATENATE_BINS.out)

    ch_bwa_mem_in = clean_fqs.combine( BWA_INDEX.out.genome_index.toList() )
    BWA_MEM(ch_bwa_mem_in)

    COVERM(BWA_MEM.out.aligned_bam, params)
    CLEAN_BAM(COVERM.out.mags_profile.join(BWA_MEM.out.aligned_bam).join(BWA_MEM.out.flag_file))

    ch_merge_coverm_in = COVERM.out.mags_profile
        .collect( flat: false, sort: { a, b -> a[0] <=> b[0] } )
        .transpose()
        .toList()
    MERGE_COVERM(ch_merge_coverm_in)


    emit:
    MERGE_COVERM.out.merged_mags_profile
}