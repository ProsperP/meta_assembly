process CLEAN_BAM {
    tag sample
    label 'process_single'

    storeDir 'cleaned_file_flags/align_MAGs_bams/'

    input:
    tuple val(sample), val(profile), val(bam), path(flag_file)

    output:
    path(flag_file)

    script:
    """
    if [ -f ${bam} ]; then
        rm ${bam}
    fi
    """
}