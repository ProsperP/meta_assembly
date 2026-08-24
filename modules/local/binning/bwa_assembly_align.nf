process bwa_mem {
    tag sample
    label 'process_medium'

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs), path(assembly)
    //val params

    output:
    tuple val(sample), eval("echo \$(pwd)/${sample}/${sample}.sorted.bam"), emit: bam
    tuple val(sample), path("${sample}/${sample}.sorted.bam.done"), emit: flag_file
    path "versions.yml", emit: versions

    script:
    """
    mkdir ${sample}/
    bwa index ${assembly}

    bwa mem -v 1 -t ${task.cpus} \\
        ${assembly} \\
        ${clean_fqs[0]} ${clean_fqs[1]} \\
        -o ${sample}/${sample}.sam

    rm ${assembly}.{amb,ann,bwt,pac,sa}

    samtools sort -T ${sample}/samtools.tmp \\
        -@ ${task.cpus} -O BAM \\
        -o ${sample}/${sample}.sorted.bam \\
        ${sample}/${sample}.sam \\
        && rm ${sample}/${sample}.sam

    touch ${sample}/${sample}.sorted.bam.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$( bwa 2>&1 | grep -oP 'Version: \\K[\\d.]+-r[\\d.]+' )
        samtools: \$( samtools 2>&1 | grep -oP 'Version: \\K[\\d.]+' )
    END_VERSIONS
    """

    stub:
    """
    mkdir ${sample}/
    touch ${sample}/${sample}.sorted.bam
    touch ${sample}/${sample}.sorted.bam.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$( bwa 2>&1 | grep -oP 'Version: \\K[\\d.]+-r[\\d.]+' )
        samtools: \$( samtools 2>&1 | grep -oP 'Version: \\K[\\d.]+' )
    END_VERSIONS
    """
}