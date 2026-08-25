process GTDBTK_CLASSIFY_WF {
    label 'process_high_memory'

    conda "/work/home/prosperp/src/miniforge3/envs/gtdbtk-2.7.2"

    input:
    path bin_files

    output:
    path "identify/gtdbtk.*.tsv", emit: identify
    path "align/gtdbtk.bac120.filtered.tsv", emit: bac120_filtered
    path "align/gtdbtk.bac120.msa.fasta.gz", emit: msa
    path "align/gtdbtk.bac120.user_msa.fasta.gz", emit: user_msa
    path "classify/gtdbtk.*.summary.tsv", emit: taxonomy

    script:
    """
    export TMPDIR=\$(mktemp -d ~/tmp/gtdbtk.XXXXXX)

    gtdbtk classify_wf \\
        --genome_dir ./ \\
        --extension fa \\
        --cpus ${task.cpus} \\
        --out_dir ./

    if [ -d \${TMPDIR} ]; then
        rm -rf \${TMPDIR}
    fi
    """

    stub:
    """
    mkdir align classify identify

    touch identify/gtdbtk.{ar53,bac120}.markers_summary.tsv
    touch identify/gtdbtk.failed_genomes.tsv
    touch identify/gtdbtk.translation_table_summary.tsv

    touch align/gtdbtk.bac120.filtered.tsv align/gtdbtk.bac120.msa.fasta.gz
    touch align/gtdbtk.bac120.user_msa.fasta.gz
    touch classify/gtdbtk.{ar53,bac120}.summary.tsv
    """
}