process COVERM {
    tag sample
    label 'process_single'

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), val(bam)
    val params

    output:
    tuple val(sample), path("${sample}_MAGs_coverm.tsv"), emit: mags_profile
    path "versions.yml", emit: versions

    script:
    """
    coverm genome \\
        --bam-files ${bam} \\
        --separator '-' \\
        --min-read-percent-identity 95 \\
        --min-read-aligned-percent 75 \\
        --min-covered-fraction 10 \\
        -t ${task.cpus} \\
        -m mean relative_abundance tpm covered_fraction \\
        -o ${sample}_MAGs_coverm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}_MAGs_coverm.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(coverm --version)
    END_VERSIONS
    """
}


process MERGE_COVERM {
    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(samples), path(profiles)

    output:
    path "all_MAGs_coverm.tsv", emit: merged_mags_profile

    script:
    """
    merge_coverm.R all_MAGs_coverm.tsv ${profiles}
    """

    stub:
    """
    touch all_MAGs_coverm.tsv
    """
}