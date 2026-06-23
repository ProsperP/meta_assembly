process checkm {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/checkm_env"

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}/checkm"), emit: checkm_out
    tuple val(sample), path("${sample}/checkm/${sample}.tsv"), emit: checkm_tsv

    script:
    """
    checkm lineage_wf \\
        -t ${task.cpus} \\
        --pplacer_threads ${task.cpus} \\
        -x ${assembly} \\

    """

}
