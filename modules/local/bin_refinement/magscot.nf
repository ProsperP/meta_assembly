process magscot {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/magscot_env"

    input:
    tuple val(sample), val(comb_name), path(bin_folders)

    output:
    tuple val(sample), val(comb_name), path(bin_folders)

    when:

}
