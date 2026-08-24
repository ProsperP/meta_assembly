process CLEAN_BAM {
    tag sample
    label 'process_single'

    storeDir 'cleaned_file_flags/align_assembly_bams/'

    input:
    tuple val(sample), val(binner), path(bin_folders), val(bam), path(flag_file)

    output:
    path "${flag_file}"

    script:
    """
    if [ -f ${bam} ]; then
        rm ${bam}
    fi
    """
}