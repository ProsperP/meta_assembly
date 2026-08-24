process QUAST {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(final_assembly)

    output:
    tuple val(sample), path("${sample}/quast/assembly_report.{html,pdf}"), emit: assembly_report
    //tuple val(sample), path("${sample}/quast/quast.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    quast -t ${task.cpus} \\
        --output-dir ${sample}/quast \\
        --min-contig 500 \\
        --silent \\
        ${final_assembly} \\
        && mv ${sample}/quast/report.html ${sample}/quast/assembly_report.html \\
        && mv ${sample}/quast/report.pdf ${sample}/quast/assembly_report.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        QUAST: \$( quast -v )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/quast
    touch ${sample}/quast/assembly_report.{html,pdf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        QUAST: \$( quast -v )
    END_VERSIONS
    """
}