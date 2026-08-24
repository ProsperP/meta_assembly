process COMEBIN {
    tag sample
    label 'process_high'

    conda "/work/home/prosperp/src/miniforge3/envs/comebin_env"

    input:
    tuple val(sample), path(assembly), val(bam)
    //val params

    output:
    tuple val(sample), val('COMEBin'), path("${sample}/comebin/comebin_bins/"), emit: bin_folder
    tuple val(sample), path("${sample}/comebin/comebin_bins/*.fa.gz"), emit: bins
    tuple val(sample), path("${sample}/comebin/comebin_res.tsv"),   emit: tsv
    tuple val(sample), path("${sample}/comebin/embeddings.tsv"),    emit: embeddings
    tuple val(sample), path("${sample}/comebin/covembeddings.tsv"), emit: covembeddings
    tuple val(sample), path("${sample}/comebin/comebin.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    mkdir ${sample}/
    pigz -d -c ${assembly} > ${sample}/${sample}_assembly.fa

    run_comebin.sh \\
        -t ${task.cpus} \\
        -a ${sample}/${sample}_assembly.fa \\
        -p \$(dirname ${bam}) \\
        -o ${sample}/ \\
        && rm ${sample}/${sample}_assembly.fa

    mv ${sample}/comebin_res ${sample}/comebin
    mv ${sample}/comebin/comebin_res_bins ${sample}/comebin/comebin_bins

    find ${sample}/comebin/comebin_bins/*.fa -exec pigz {} \\;

    rm -rf ${sample}/data_augmentation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | grep -oP 'version: \\K[\\d.]+')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/comebin/comebin_bins

    echo '' | gzip > ${sample}/comebin/comebin_bins/1.fa.gz
    echo '' | gzip > ${sample}/comebin/comebin_bins/2.fa.gz

    touch ${sample}/comebin/comebin_res.tsv
    touch ${sample}/comebin/{embeddings,covembeddings}.tsv
    touch ${sample}/comebin/comebin.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | grep -oP 'version: \\K[\\d.]+')
    END_VERSIONS
    """
}