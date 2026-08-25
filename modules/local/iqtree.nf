process IQTREE3 {
    label 'process_proper'

    conda "/work/home/prosperp/src/miniforge3/envs/gtdbtk-2.7.2"

    input:
    path user_msa

    output:
    path "tree/gtdbtk.bac120.*"

    script:
    """
    mkdir tree

    n_seqs=\$( zcat ${user_msa} | grep '>' | wc -l )
    if [ \${n_seqs} -ge 4 ]; then
        iqtree3 \\
            -s ${user_msa} \\
            -m MFP --ufboot 1000 -alrt 1000 \\
            -T ${task.cpus} \\
            --prefix tree/gtdbtk.bac120
    else
        iqtree3 \\
            -s ${user_msa} \\
            -m MFP -alrt 1000 \\
            -T ${task.cpus} \\
            --prefix tree/gtdbtk.bac120
    fi
    """

    stub:
    """
    mkdir tree
    touch tree/gtdbtk.bac120.{bionj,contree,log,model.gz,treefile}
    touch tree/gtdbtk.bac120.{ckp.gz,iqtree,mldist,splits.nex}
    """
}