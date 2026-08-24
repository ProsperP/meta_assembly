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
    path "versions.yml", emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaSPAdes: \$( metaspades.py -v | grep -oP 'v[\\d.]+' )
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaSPAdes: \$( metaspades.py -v | grep -oP 'v[\\d.]+' )
    END_VERSIONS
    """
}


process megahit {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(clean_fqs)
    val min_len

    output:
    tuple val(sample), path("${sample}/${sample}.final_assembly.fa.gz"), emit: final_assembly
    tuple val(sample), path("${sample}/log"), emit: log
    path "versions.yml", emit: versions

    script:
    def mem = task.memory.toBytes()
    """
    mkdir -p ${sample}_megahit.tmp

    megahit \\
        -1 ${clean_fqs[0]} -2 ${clean_fqs[1]} \\
        -o ${sample}/ \\
        --tmp-dir ${sample}_megahit.tmp \\
        -t ${task.cpus} -m ${mem} \\
        --continue

    if [ -d ${sample}_megahit.tmp ]; then
        rm -rf ${sample}_megahit.tmp
    fi

    fix_megahit_contig_naming.py \\
        ${sample}/final.contigs.fa ${min_len} \\
        | gzip > ${sample}/${sample}.final_assembly.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$( megahit -v )
    END_VERSIONS
    """

    stub:
    """
    mkdir ${sample}/
    touch ${sample}/${sample}.final_assembly.fa.gz ${sample}/log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$( megahit -v )
    END_VERSIONS
    """
}


process filterAssembly {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(scaffolds)
    val min_len

    output:
    tuple val(sample), path("${sample}/${sample}.final_assembly.fa.gz"), emit: final_assembly

    script:
    """
    mkdir ${sample}/
    seqkit seq --min-len ${min_len} ${scaffolds} \\
        | gzip > ${sample}/${sample}.final_assembly.fa.gz
    """

    stub:
    """
    mkdir ${sample}/
    touch ${sample}/${sample}.final_assembly.fa.gz
    """
}