process DREP {
    label 'process_high'

    conda "/work/home/prosperp/src/miniforge3/envs/drep_env"

    input:
    path bin_files
    val min_ani

    output:
    path "dereplicated_genomes/*.fa"

    script:
    """
    OUTDIR=./

    mkdir -p \${OUTDIR}

    dRep dereplicate \${OUTDIR} \\
        -sa ${min_ani} \\
        -comp 50 -con 10 \\
        -p ${task.cpus} \\
        -g ${bin_files}


    """

    stub:
    def outfiles = bin_files.collect { it -> "dereplicated_genomes/${it.name}" }.join(' ')
    """
    mkdir dereplicated_genomes
    touch ${outfiles}
    """
}


process RENAME_BINS {
    label 'process_single'

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    path bin_files

    output:
    path "bins_mags_id.txt", emit: id_convert
    path "MAGs/*.fa",        emit: mags

    script:
    """
    rename_mags.py -g ./ -o ./
    mv MAGs_newid MAGs
    """

    stub:
    """
    rename_mags.py -g ./ -o ./
    mv MAGs_newid MAGs
    """
}