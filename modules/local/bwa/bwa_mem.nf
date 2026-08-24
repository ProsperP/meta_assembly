process BWA_MEM {
    tag sample
    label 'process_medium'
    label 'process_medium_memory'

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs), path(index)

    output:
    tuple val(sample), eval("echo \$(pwd)/${sample}.bam"), emit: aligned_bam
    tuple val(sample), path("${sample}.bam.done"), emit: flag_file
    path "versions.yml", emit: versions

    script:
    def index_name = index[0].baseName
    """
    bwa mem -v 1 -t ${task.cpus} \\
        ${index_name} \\
        ${clean_fqs} \\
        -o ${sample}.sam

    samtools sort -T ${sample}_samtools.tmp \\
        -@ ${task.cpus} -O BAM \\
        -o ${sample}.bam \\
        ${sample}.sam \\
        && rm ${sample}.sam

    touch ${sample}.bam.done

    if [ -d ${sample}_samtools.tmp ]; then
        rm -rf ${sample}_samtools.tmp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$( bwa 2>&1 | grep -oP 'Version: \\K[\\d.]+-r[\\d.]+' )
        samtools: \$( samtools 2>&1 | grep -oP 'Version: \\K[\\d.]+' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.bam
    touch ${sample}.bam.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$( bwa 2>&1 | grep -oP 'Version: \\K[\\d.]+-r[\\d.]+' )
        samtools: \$( samtools 2>&1 | grep -oP 'Version: \\K[\\d.]+' )
    END_VERSIONS
    """
}