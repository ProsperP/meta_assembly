process summarize_checkm {
    tag "${sample}-${comb_name}"

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), val(comb_name), path(checkm2_tsv)

    output:
    tuple val(sample), val(comb_name), path("${sample}/${sample}_${comb_name}.stats"), emit: checkm2_stats

    script:
    """
    mkdir -p ${sample}/
    summarize_checkm.py -i ${checkm2_tsv} -l ${comb_name} \\
        | (read -r; printf "%s\n" "\$REPLY"; sort) \\
        > ${sample}/${sample}_${comb_name}.stats
    """

    stub:
    """
    mkdir -p ${sample}/
    touch ${sample}/${sample}_${comb_name}.stats
    echo summarize_checkm.py -i ${checkm2_tsv} -l ${comb_name}
    """
}
