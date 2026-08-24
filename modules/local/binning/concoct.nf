process CONCOCT {
    tag sample
    label 'process_medium'

    conda "/work/home/prosperp/src/miniforge3/envs/concoct_env"

    input:
    tuple val(sample), path(assembly), val(bam)
    val params

    output:
    tuple val(sample), val('CONCOCT'), path("${sample}/concoct/concoct_bins/"), emit: bin_folder
    tuple val(sample), path("${sample}/concoct/concoct_bins/*.fa.gz"), emit: bins
    tuple val(sample), path("${sample}/concoct/concoct_depth.tsv"), emit: concoct_depth
    tuple val(sample), path("${sample}/concoct/log.txt"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    outdir=${sample}/concoct
    mkdir -p \${outdir}

    pigz -d -c ${assembly} > ${sample}/${sample}_assembly.fa

    samtools index -@ ${task.cpus} -b ${bam}

    cut_up_fasta.py \\
        ${sample}/${sample}_assembly.fa \\
        --chunk_size 10000 \\
        --merge_last --overlap_size 0 \\
        --bedfile \${outdir}/assembly_10K.bed \\
        > \${outdir}/assembly_10K.fa

    concoct_coverage_table.py \\
        \${outdir}/assembly_10K.bed \\
        ${bam} \\
        > \${outdir}/concoct_depth.tsv

    # Starting binning with CONCOCT...
    concoct \\
        --threads ${task.cpus} \\
        -l ${params.min_len} \\
        --coverage_file \${outdir}/concoct_depth.tsv \\
        --composition_file \${outdir}/assembly_10K.fa \\
        --basename \${outdir}/

    merge_cutup_clustering.py \\
        \${outdir}/clustering_gt${params.min_len}.csv \\
        > \${outdir}/clustering_gt${params.min_len}_merged.csv

    # splitting contigs into bins
    mkdir \${outdir}/concoct_bins
    extract_fasta_bins.py \\
        ${sample}/${sample}_assembly.fa \\
        \${outdir}/clustering_gt${params.min_len}_merged.csv \\
        --output_path \${outdir}/concoct_bins

    find \${outdir}/concoct_bins/*.fa -exec pigz {} \\;

    # Cleanup
    rm ${bam}.bai
    rm ${sample}/${sample}_assembly.fa
    rm \${outdir}/assembly_10K* \${outdir}/*.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONCOCT: \$(concoct -v)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/concoct/concoct_bins

    echo '' | gzip > ${sample}/concoct/concoct_bins/{1,2,3,4,5}.fa.gz
    touch ${sample}/concoct/log.txt ${sample}/concoct/concoct_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONCOCT: \$(concoct -v)
    END_VERSIONS
    """
}