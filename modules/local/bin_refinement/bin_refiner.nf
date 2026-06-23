process bin_refiner {
    tag "${sample}-${comb_name}"

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), val(comb_name), path(bin_folders)

    output:
    tuple val(sample), val(comb_name), path("${sample}/${comb_name}_bins"), emit: bins_comb

    when:
    bin_folders.size() > 1

    script:
    def binset3_opt = bin_folders.size() == 3 ? "-3 ${sample}/${bin_folders[2]}" : ''
    """
    mkdir -p ${sample}/
    for dir in ${bin_folders}; do
        if [[ -d \${dir} ]]; then
            mkdir -p ${sample}/\${dir}
            parallel -j ${task.cpus} '
                SIZE=\$(pigz -l {} | awk "NR==2 {print \\\$2}")
                if (( \$SIZE > 50000 )) && (( \$SIZE < 20000000 )); then
                    pigz -d -c {} | sed '/^>/ s/=/_/g' > ${sample}/\${dir}/{/.}
                else
                    echo "Skipping {} since its bin size if not between 50kb and 20Mb"
                fi
                ' ::: \$(ls \${dir}/*.fa.gz)
        else
            error "\${dir} is not a valid directory."
        fi
    done

    binning_refiner.py \\
        -1 ${sample}/${bin_folders[0]} \\
        -2 ${sample}/${bin_folders[1]} \\
        ${binset3_opt} \\
        -o ${sample}/Refined_${comb_name}

    mv ${sample}/Refined_${comb_name}/Refined_bins ${sample}/${comb_name}_bins

    # Clean up
    rm -rf ${sample}/Refined_${comb_name}
    """

    stub:
    def binset3_opt = bin_folders.size() == 3 ? "-3 ${sample}/${bin_folders[2]}" : ''
    """
    for dir in ${bin_folders}; do
        if [[ -d \${dir} ]]; then
            mkdir -p ${sample}/\${dir}
            parallel -j 1 '
                SIZE=\$(pigz -l {} | awk "NR==2 {print \\\$2}")
                if (( \${SIZE} > 0 )); then
                    echo \$(pigz -l {} | awk "NR==2 {print \\\$2}")
                fi
                ' ::: \$(ls \${dir}/*)
        fi
    done
    echo binning_refiner.py \\
        -1 ${sample}/${bin_folders[0]} \\
        -2 ${sample}/${bin_folders[1]} \\
        ${binset3_opt} \\
        -o ${sample}/Refined_${comb_name}
    mkdir -p ${sample}/${comb_name}_bins
    """
}
