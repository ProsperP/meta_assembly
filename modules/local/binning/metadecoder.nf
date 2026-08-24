process METADECODER {
    tag sample
    label 'process_proper'

    conda "/work/home/prosperp/src/miniforge3/envs/metadecoder_env"

    input:
    tuple val(sample), path(assembly), val(bam)
    val params

    output:
    tuple val(sample), val('MetaDecoder'), path("${sample}/metadecoder/metadecoder_bins/"), emit: bin_folder
    tuple val(sample), path("${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder.*.fa.gz"), emit: bins
    path "versions.yml", emit: versions

    script:
    """
    mkdir ${sample}
    metadecoder coverage \\
        --bam ${bam} \\
        --output ${sample}/${sample}.metadecoder.coverage

    metadecoder seed \\
        --threads ${task.cpus} \\
        --fasta ${assembly} \\
        --output ${sample}/${sample}.metadecoder.seed

    mkdir -p ${sample}/metadecoder/metadecoder_bins
    metadecoder cluster \\
        ${params.ext_args} \\
        --fasta ${assembly} \\
        --coverage ${sample}/${sample}.metadecoder.coverage \\
        --seed ${sample}/${sample}.metadecoder.seed \\
        --output ${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder

    rename .fasta .fa ${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder.*.fasta

    rm ${sample}/${sample}.metadecoder.{coverage,seed} ${bam}.index

    find ${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder.*.fa -exec pigz {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaDecoder: \$(metadecoder -v)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/metadecoder/metadecoder_bins
    touch ${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder.{1,2,3,4,5}.fa
    gzip ${sample}/metadecoder/metadecoder_bins/${sample}.metadecoder.*.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaDecoder: \$(metadecoder -v)
    END_VERSIONS
    """
}