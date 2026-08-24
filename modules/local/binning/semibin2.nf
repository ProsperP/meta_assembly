process SEMIBIN2 {
    tag sample
    label 'process_low'

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(assembly), val(bam)
    val params

    output:
    tuple val(sample), path("${sample}/semibin2/*.csv"), emit: csv
    tuple val(sample), path("${sample}/semibin2/*.tsv"), emit: tsv
    tuple val(sample), path("${sample}/semibin2/*.h5"),  emit: model, optional: true
    tuple val(sample), path("${sample}/semibin2/semibin2_bins/*.fa.gz"),      emit: bins
    tuple val(sample), val('SemiBin2'), path("${sample}/semibin2/semibin2_bins/"),      emit: bin_folder
    tuple val(sample), path("${sample}/semibin2/${sample}_SemiBinRun.log"), emit: log
    path "versions.yml", emit: versions


    script:
    """
    mkdir -p ${sample}/semibin2

    SemiBin2 single_easy_bin \\
        --threads ${task.cpus} \\
        --self-supervised \\
        ${params.ext_args} \\
        --input-fasta ${assembly} \\
        --input-bam ${bam} \\
        --output ${sample}/semibin2 \\
        --tmpdir ${sample}/semibin2.tmp

    mv ${sample}/semibin2/output_bins ${sample}/semibin2/semibin2_bins
    rmdir ${sample}/semibin2.tmp

    mv ${sample}/semibin2/SemiBinRun.log ${sample}/semibin2/${sample}_SemiBinRun.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$(SemiBin2 --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/semibin2/semibin2_bins

    touch ${sample}/semibin2/{contig_bins,recluster_bins_info}.tsv
    touch ${sample}/semibin2/{data,data_split}.csv
    touch ${sample}/semibin2/${sample}_SemiBinRun.log

    echo '' | gzip > ${sample}/semibin2/semibin2_bins/1.fa.gz
    echo '' | gzip > ${sample}/semibin2/semibin2_bins/2.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SemiBin2: \$(SemiBin2 --version)
    END_VERSIONS
    """
}