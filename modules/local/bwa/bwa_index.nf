process BWA_INDEX {

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    path genome

    output:
    path "${genome}.*", emit: genome_index

    script:
    """
    bwa index ${genome} -p ./${genome}
    """
}