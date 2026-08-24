process collectFlagFiles {
    label 'process_single'

    storeDir('temp_file_flag/')

    input:
    tuple val(samples), path(flag_files)
    val file_name

    output:
    path "${file_name}"

    script:
    """
    > ${file_name}
    for flag_file in ${flag_files}; do
        echo \$(readlink -f \${flag_file}) >> ${file_name}
    done
    """

    stub:
    """
    touch ${file_name}
    """
}