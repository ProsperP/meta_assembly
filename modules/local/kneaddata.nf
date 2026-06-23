process kneaddata {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/biobakery3"

    input:
    tuple val(sample), path(input_fqs)
    val params

    output:
    tuple val(sample), path("${sample}/${sample}_paired_{1,2}.fastq", arity: '2'), emit: clean_fqs
    tuple val(sample), path("${sample}/${sample}_kneaddata.log"), emit: log

    when:
    !params.skip

    script:
    trim_mem = "6g"
    """
    mkdir ${sample}/
    export _JAVA_OPTIONS="-Xmx${trim_mem}"
    kneaddata --input1 ${input_fqs[0]} --input2 ${input_fqs[1]} \\
        --reference-db ${params.bowtie2_db} \\
        --output ${sample}/ --output-prefix ${sample} \\
        --max-memory ${trim_mem} \\
        --trimmomatic-options "ILLUMINACLIP:${params.adapter}:2:30:10 LEADING:${params.lead_qual} TRAILING:${params.trail_qual} SLIDINGWINDOW:${params.window_size}:${params.avg_quality} HEADCROP:${params.head_crop} MINLEN:${params.min_len}" \\
        --remove-intermediate-output --run-fastqc-start --run-fastqc-end \\
        --threads ${task.cpus} --log ${sample}/${sample}_kneaddata.log \\
        && gzip ${sample}/*_unmatched_?.fastq ${sample}/*_contam*fastq
    """

    stub:
    """
    mkdir ${sample}/
    kneaddata -h > ${sample}/${sample}_kneaddata.log
    ls ${input_fqs[0]} ${input_fqs[1]} \\
        && touch '${sample}/${sample}_paired_1.fastq' '${sample}/${sample}_paired_2.fastq'
    """
}
