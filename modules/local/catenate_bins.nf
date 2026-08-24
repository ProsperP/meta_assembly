process CATENATE_BINS {
    label 'process_single'

    input:
    path bin_files

    output:
    path "all_MAGs.fa"

    script:
    """
    cat \$(ls -v ./) > all_MAGs.fa
    """

    stub:
    """
    touch all_MAGs.fa
    """
}