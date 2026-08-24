include { DREP as DREP_BATCH } from '../modules/local/drep'
include { DREP as DREP_MAGS } from '../modules/local/drep'
include { DREP as DREP_SRGS } from '../modules/local/drep'
include { RENAME_BINS } from '../modules/local/drep'


workflow DEREPLICATION {
    take:
    ch_all_bins
    params
    drep_batch_size

    main:

    ch_all_samples = ch_all_bins
        .map { _sample, list -> list }
        .collect( flat: false )
        .branch { list ->
            split: list.size() >= drep_batch_size
                return list
            single: list.size() < drep_batch_size
                return list
        }

    ch_in_drep_batches = ch_all_samples.split
        .map { list -> list.collate( drep_batch_size ) }
        .flatMap()
        .map { list -> list.flatten() }

    DREP_BATCH(ch_in_drep_batches, params.mags_ani)

    ch_final_result = ch_all_samples.single
        .map { list -> list.flatten() }
        .mix(DREP_BATCH.out.collect())

    DREP_MAGS(ch_final_result, params.mags_ani)
    RENAME_BINS(DREP_MAGS.out)
    ch_mags = RENAME_BINS.out.mags
    DREP_SRGS(ch_mags, params.srgs_ani)

    emit:
    id_convert = RENAME_BINS.out.id_convert
    mags = ch_mags
    srgs = DREP_SRGS.out
}