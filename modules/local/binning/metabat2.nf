process metabat2 {
    tag sample

    conda "/work/home/prosperp/src/miniforge3/envs/meta_assembly"

    input:
    tuple val(sample), path(assembly), path(bam)
    val params

    output:
    tuple val(sample), path("${sample}/metabat2/*.tooShort.fa.gz"), emit: tooshort, optional: true
    tuple val(sample), path("${sample}/metabat2/*.lowDepth.fa.gz"), emit: lowdepth, optional: true
    tuple val(sample), path("${sample}/metabat2/*.unbinned.fa.gz"), emit: unbinned, optional: true
    tuple val(sample), path("${sample}/metabat2/metabat2_bins/*[!lowDepth|tooShort|unbinned].fa.gz"),  emit: bins, optional: true
    tuple val(sample), val('MataBAT2'), path("${sample}/metabat2/metabat2_bins/"),  emit: bin_folder, optional: true
    tuple val(sample), path("${sample}/metabat2/assembly_depth.txt"), emit: depth
    //tuple val(sample), path("*.tsv.gz"), emit: membership, optional: true
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p ${sample}/metabat2

    jgi_summarize_bam_contig_depths \\
        --outputDepth ${sample}/metabat2/assembly_depth.txt \\
        ${bam}

    metabat2 \\
        --numThreads ${task.cpus} \\
        --unbinned \\
        --minContig ${params.min_len} \\
        --inFile ${assembly} \\
        --abdFile ${sample}/metabat2/assembly_depth.txt \\
        --outFile ${sample}/metabat2/${sample}

    mkdir ${sample}/metabat2/metabat2_bins/
    pigz ${sample}/metabat2/${sample}.*.fa

    shopt -s extglob
    if [[ -f ${sample}/metabat2/${sample}.1.fa.gz ]]; then
        mv ${sample}/metabat2/!(*lowDepth|*tooShort|*unbinned).fa.gz ${sample}/metabat2/metabat2_bins/

        cd ${sample}/metabat2/metabat2_bins && rename '${sample}.' '' *.fa.gz && cd -
    fi
    shopt -u extglob

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 -h 2>&1 | grep -oP 'version 2:\\K[\\d.]+' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample}/metabat2/metabat2_bins/

    touch ${sample}/metabat2/assembly_depth.txt

    echo "" | gzip -c > ${sample}/metabat2/${sample}.1.fa.gz
    echo "" | gzip -c > ${sample}/metabat2/${sample}.2.fa.gz
    echo "" | gzip -c > ${sample}/metabat2/${sample}.1.tooShort.fa.gz
    echo "" | gzip -c > ${sample}/metabat2/${sample}.1.lowDepth.fa.gz
    echo "" | gzip -c > ${sample}/metabat2/${sample}.1.unbinned.fa.gz
    #echo "" | gzip -c > ${sample}/metabat2/${sample}.tsv.gz

    shopt -s extglob
    mv ${sample}/metabat2/!(*lowDepth|*tooShort|*unbinned).fa.gz ${sample}/metabat2/metabat2_bins/
    shopt -u extglob

    cd ${sample}/metabat2/metabat2_bins && rename '${sample}.' '' *.fa.gz && cd -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabat2: \$( metabat2 -h 2>&1 | grep -oP 'version 2:\\K[\\d.]+' )
    END_VERSIONS
    """
}
