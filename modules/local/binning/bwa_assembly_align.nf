process bwa_mem {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs), path(assembly)
    //val params

    output:
    tuple val(sample), path("${sample}/${sample}.sorted.bam"), emit: bam
    tuple val(sample), path("${sample}/${sample}.bwa.log"), emit: log

    script:
    """
    mkdir ${sample}/
    {
        bwa index ${assembly}

        bwa mem -v 1 -t ${task.cpus} \\
            ${assembly} \\
            ${clean_fqs[0]} ${clean_fqs[1]} \\
            -o ${sample}/${sample}.sam

        samtools sort -T ${sample}/samtools.tmp \\
            -@ ${task.cpus} -O BAM \\
            -o ${sample}/${sample}.sorted.bam \\
            ${sample}/${sample}.sam \\
            && rm ${sample}/${sample}.sam
    } 2>&1 > ${sample}/${sample}.bwa.log
    """

    stub:
    """
    mkdir ${sample}/
    touch ${sample}/${sample}.sorted.bam ${sample}/${sample}.bwa.log
    """
}
