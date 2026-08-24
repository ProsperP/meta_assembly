process identify_single_copy_marker {
    tag sample
    label 'process_low'

    conda "/work/home/prosperp/src/miniforge3/envs/magscot_env"

    input:
    tuple val(sample), path(assembly)
    val params

    output:
    tuple val(sample), path("${sample}.prodigal.faa"), emit: faa
    tuple val(sample), path("${sample}.prodigal.ffn"), emit: ffn
    tuple val(sample), path("${sample}.hmm"), emit: hmm
    path "versions.yml", emit: versions

    script:
    """
    tmpdir=tmp_workfolder
    mkdir \${tmpdir}

    zcat ${assembly} | parallel -j ${task.cpus} --block 999k --recstart '>' --pipe prodigal -p meta \\
        -a \${tmpdir}/${sample}.{#}.faa -d \${tmpdir}/${sample}.{#}.ffn \\
        -o tmpfile

    cat \${tmpdir}/${sample}.*.faa > ${sample}.prodigal.faa
    cat \${tmpdir}/${sample}.*.ffn > ${sample}.prodigal.ffn
    rm -r \${tmpdir} tmpfile

    hmmsearch \\
        -o ${sample}.hmm.tigr.out \\
        --tblout ${sample}.hmm.tigr.hit.out \\
        --noali --notextw --cut_nc \\
        --cpu ${task.cpus} \\
        ${params.tigr_db} \\
        ${sample}.prodigal.faa

    hmmsearch \\
        -o ${sample}.hmm.pfam.out \\
        --tblout ${sample}.hmm.pfam.hit.out \\
        --noali --notextw --cut_nc \\
        --cpu ${task.cpus} \\
        ${params.pfam_db} \\
        ${sample}.prodigal.faa

    cat ${sample}.hmm.tigr.hit.out | grep -v "^#" | awk '{print \$1"\\t"\$3"\\t"\$5}' > ${sample}.tigr
    cat ${sample}.hmm.pfam.hit.out | grep -v "^#" | awk '{print \$1"\\t"\$4"\\t"\$5}' > ${sample}.pfam
    cat ${sample}.tigr ${sample}.pfam > ${sample}.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Prodigal: \$(prodigal -v 2>&1 | grep -oP 'Prodigal V\\K[\\d.]+' )
        hmmsearch: \$( hmmsearch -h | grep -oP '^# HMMER \\K[\\d.]+' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.prodigal.faa ${sample}.prodigal.ffn
    touch ${sample}.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Prodigal: \$(prodigal -v 2>&1 | grep -oP 'Prodigal V\\K[\\d.]+' )
        hmmsearch: \$( hmmsearch -h | grep -oP '^# HMMER \\K[\\d.]+' )
    END_VERSIONS
    """
}