process magscot {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/magscot_env"

    input:
    tuple val(sample), path(assembly), path(bin_folders), path(hmm)
    val params

    output:
    tuple val(sample), path("MAGScoT.refined.contig_to_bin.out")
    tuple val(sample), path("${sample}/refined_bins/*.fa")


    script:
    """
    Rscript ${params.MAGScoT_folder}/MAGScoT.R \\
        -i ${sample}/contigs_to_bins.tsv \\
        -o ${sample}/MAGScoT \\
        --max_cont 0.1 \\
        --hmm ${sample}/${sample}.hmm

    mkdir ${sample}/refined_bins
    extract_refined_bins.py \\
        ${assembly} \\
        ${sample}/MAGScoT.refined.contig_to_bin.out \\
        --output_path ${sample}/refined_bins
    """

    stub:
    """
    touch ${sample}/MAGScoT.refined.contig_to_bin.out
    mkdir ${sample}/refined_bins
    touch ${sample}/refined_bins/${sample}.fa
    """

}
