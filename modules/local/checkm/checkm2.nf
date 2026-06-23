process checkm2 {
    tag "${sample}-${comb_name}"

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), val(comb_name), path(bin_set)
    val params

    output:
    tuple val(sample), path("${sample}/checkm2"), emit: checkm2_out
    tuple val(sample), val(comb_name), path("${sample}/checkm2/${sample}_${comb_name}_checkm2_report.tsv"), emit: checkm2_tsv
    path "versions.yml", emit: versions

    script:
    """
    checkm2 predict \\
        --threads ${task.cpus} \\
        --extension fa.gz \\
        --input ${bin_set} \\
        --output-directory ${sample}/checkm2 \\
        --tmpdir ${sample}/checkm2/${sample}.tmp \\
        --database_path ${params.db}

    cp ${sample}/checkm2/quality_report.tsv ${sample}/checkm2/${sample}_${comb_name}_checkm2_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkM2: \$(checkm2 --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/checkm2
    touch ${sample}/checkm2/${sample}_${comb_name}_checkm2_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkM2: \$(checkm2 --version)
    END_VERSIONS
    """
}
