process format_contigs2bin {
    tag "${sample}-${binner}"
    label 'process_single'

    conda "/work/home/prosperp/src/miniforge3/envs/magscot_env"

    input:
    tuple val(sample), val(binner), path(bin_folders)

    output:
    tuple val(sample), val(binner), path("${binner}.contigs_to_bin.tsv")

    script:
    """
    contigs_to_bin_tsv.py \\
        -p ${bin_folders} \\
        -o ${binner}.contigs_to_bin.tsv \\
        -l ${binner}
    """

    stub:
    """
    touch ${binner}.contigs_to_bin.tsv
    """
}


process magscot {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/magscot_env"

    input:
    tuple val(sample), path(assembly), path(contigs2bin_tsvs), path(hmm)
    val params

    output:
    tuple val(sample), path("${sample}/MAGScoT.*.out"), emit: magscot_out
    tuple val(sample), path("${sample}/refined_bins/${sample}_MAGScoT*.fa"), emit: refined_bins

    script:
    """
    mkdir ${sample}/
    cat ${contigs2bin_tsvs} > ${sample}/contigs_to_bin.tsv

    Rscript ${params.magscot_dir}/MAGScoT.R \\
        -i ${sample}/contigs_to_bin.tsv \\
        -o ${sample}/MAGScoT \\
        --max_cont ${params.contamination} \\
        --hmm ${hmm}

    mkdir ${sample}/refined_bins
    extract_refined_bins.py \\
        ${assembly} \\
        ${sample}/MAGScoT.refined.contig_to_bin.out \\
        --output_path ${sample}/refined_bins
    rename 'MAGScoT' '${sample}_MAGScoT' ${sample}/refined_bins/*.fa
    """

    stub:
    """
    mkdir ${sample}
    cat ${contigs2bin_tsvs} > ${sample}/contigs_to_bin.tsv

    echo Rscript ${params.magscot_dir}/MAGScoT.R \\
        -i ${sample}/contigs_to_bin.tsv \\
        -o ${sample}/MAGScoT \\
        --max_cont ${params.contamination} \\
        --hmm ${sample}/${sample}.hmm

    touch ${sample}/MAGScoT.refined.contig_to_bin.out
    touch ${sample}/MAGScoT.{refined,scores}.out

    mkdir ${sample}/refined_bins
    touch ${sample}/refined_bins/MAGScoT_cleanbin_00000{1..3}.fa
    rename 'MAGScoT' '${sample}_MAGScoT' ${sample}/refined_bins/*.fa
    """

}