process spades {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs)
    val params

    output:
    tuple val(sample), path("${sample}/*.contigs.fa.gz"),   emit: contigs
    tuple val(sample), path("${sample}/*.scaffolds.fa.gz"), emit: scaffolds
    tuple val(sample), path("${sample}/*.gene_clusters.fa.gz"), optional: true, emit: gene_clusters
    tuple val(sample), path("${sample}/*.transcripts.fa.gz"),   optional: true, emit: transcripts
    tuple val(sample), path("${sample}/*.assembly.gfa.gz"),     optional: true, emit: gfa
    tuple val(sample), path("${sample}/*.warnings.log"),        optional: true, emit: warnings
    tuple val(sample), path("${sample}/*.spades.log"), emit: log

    when:
    !params.skip

    script:
    def mem = task.memory.toGiga()
    """
    if [[ -s ${sample}/${sample}.spades.log ]]; then
        metaspades.py \\
            -o ${sample}/ \\
            --tmp-dir ${sample}/metaspades.tmp \\
            --restart-from last \\
            -t ${task.cpus} -m ${mem}
    else
        metaspades.py \\
            -o ${sample} \\
            --tmp-dir ${sample}/metaspades.tmp \\
            -1 ${clean_fqs[0]} -2 ${clean_fqs[1]} \\
            -t ${task.cpus} -m ${mem}
    fi

    rm -r ${sample}/metaspades.tmp

    cd ${sample}/
    mv spades.log ${sample}.spades.log
    for prefix in {contigs,scaffolds,gene_clusters,transcripts}; do
        if [ -f \${prefix}.fasta ]; then
            mv \${prefix}.fasta ${sample}.\${prefix}.fa \\
                && gzip -n ${sample}.\${prefix}.fa
        fi
    done

    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${sample}.assembly.gfa \\
            && gzip -n ${sample}.assembly.gfa
    fi

    if [ -f warnings.log ]; then
        mv warnings.log ${sample}.warnings.log
    fi
    cd ..
    """

    stub:
    """
    mkdir ${sample}/
    cd ${sample}/
    for prefix in {contigs,scaffolds,gene_clusters,transcripts}; do
        echo "" | gzip > ${sample}.\${prefix}.fa.gz
    done
    touch ${sample}.assembly.gfa ${sample}.warnings.log ${sample}.spades.log
    cd ..
    """
}


process filter_assembly {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(scaffolds)
    val params

    output:
    tuple val(sample), path("${sample}/${sample}.final_assembly.fa.gz"), emit: final_assembly
    tuple val(sample), path("${sample}/${sample}.filter.log"), emit: log

    script:
    """
    mkdir ${sample}/
    { seqkit seq --min-len ${params.min_len} ${scaffolds} \\
        | gzip > ${sample}/${sample}.final_assembly.fa.gz
        } 2> ${sample}/${sample}.filter.log
    """

    stub:
    """
    mkdir ${sample}/
    touch ${sample}/${sample}.final_assembly.fa.gz
    touch ${sample}/${sample}.filter.log
    """
}

process quast {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(final_assembly)

    output:
    tuple val(sample), path("${sample}/quast/assembly_report.{html,pdf}"), emit: assembly_report
    //tuple val(sample), path("${sample}/quast/quast.log"), emit: log

    script:
    """
    quast -t ${task.cpus} \\
        --output-dir ${sample}/quast \\
        --min-contig 500 \\
        --silent \\
        ${final_assembly} \\
        && mv ${sample}/quast/report.html ${sample}/quast/assembly_report.html \\
        && mv ${sample}/quast/report.pdf ${sample}/quast/assembly_report.pdf
    """

    stub:
    """
    mkdir -p ${sample}/quast
    touch ${sample}/quast/assembly_report.{html,pdf}
    """
}
