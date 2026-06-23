process concoct {
    tag sample
    label 'process_high'

    conda "/work/home/prosperp/src/miniforge3/envs/concoct_env"

    input:
    tuple val(sample), path(assembly), path(bam)
    //val params

    output:
    tuple val(sample), val('CONCOCT'), path("${sample}/concoct/concoct_bins/"), emit: bin_folder
    tuple val(sample), path("${sample}/concoct/concoct_bins/*.fa.gz"), emit: bins
    tuple val(sample), path("${sample}/concoct/comebin_res.tsv"),   emit: tsv
    tuple val(sample), path("${sample}/concoct/embeddings.tsv"),    emit: embeddings
    tuple val(sample), path("${sample}/concoct/covembeddings.tsv"), emit: covembeddings
    tuple val(sample), path("${sample}/concoct/concoct.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    mkdir ${sample}/
    pigz -d -c ${assembly} > ${sample}/${sample}_assembly.fa

    # indexing .bam alignment files...
    samtools index -@ ${task.cpus} -b ${bam}

    # cutting up contigs into 10kb fragments for CONCOCT...
    cut_up_fasta.py \\
        ${sample}/${sample}_assembly.fa \\
        -c 10000 --merge_last -o 0 \\
        -b assembly_10K.bed \\
        > ${sample}/${sample}_assembly_10K.fa

    # estimating contig fragment coverage...
    concoct_coverage_table.py \\
        assembly_10K.bed \\
        ${out}/work_files/*.bam \\
        > ${out}/work_files/concoct_depth.txt

    # Starting binning with CONCOCT...
    concoct -l $len -t $threads \\
        --coverage_file ${out}/work_files/concoct_depth.txt \\
        --composition_file ${out}/work_files/assembly_10K.fa \\
        -b ${out}/work_files/concoct_out

    # merging 10kb fragments back into contigs
    merge_cutup_clustering.py \\
        ${out}/work_files/concoct_out/clustering_gt${len}.csv \\
        > ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv

    # splitting contigs into bins
    split_concoct_bins.py \\
        ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv \\
        ${out}/work_files/assembly.fa \\
        ${out}/concoct_bins

    find ${sample}/comebin/comebin_bins/*.fa -exec pigz {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONCOCT: \$(concoct -v)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/concoct/concoct_bins

    echo '' | gzip > ${sample}/concoct/concoct_bins/1.fa.gz
    echo '' | gzip > ${sample}/concoct/concoct_bins/2.fa.gz

    touch ${sample}/comebin/comebin_res.tsv
    touch ${sample}/comebin/{embeddings,covembeddings}.tsv
    touch ${sample}/concoct/concoct.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CONCOCT: \$(concoct -v)
    END_VERSIONS
    """
}
