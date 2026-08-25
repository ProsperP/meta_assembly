process FASTQ_SORT {
    tag sample
    label "process_medium"

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs)

    output:
    tuple val(sample), path("${sample}/*.gz")

    script:
    """
    temp_dir=\$(mktemp -d ./tmpdir.XXXXXX)
    mkdir ${sample}
    for fq in ${clean_fqs}; do
        pigz -d -c \${fq} | fastq-sort --temporary-directory=\${temp_dir} \\
            | pigz --processes ${task.cpus} -c > ${sample}/\${fq}
    done

    if [ -d \${temp_dir} ]; then
        rm -rf \${temp_dir}
    fi
    """

    stub:
    """
    mkdir ${sample}
    for fq in ${clean_fqs}; do
        touch ${sample}/\${fq}
    done
    """
}