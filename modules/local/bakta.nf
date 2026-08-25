process BAKTA {
    label 'process_medium'

    conda "/work/home/prosperp/src/miniforge3/envs/bakta"

    input:
    path bins
    val params

    output:
    path "MAG_*/MAG_*.{embl,ffn,gbff,gff3,faa,fna,hypotheticals.faa,hypotheticals.tsv,inference.tsv,png,svg,tsv,txt,log}"

    script:
    """
    for bin in ${bins}; do
        bin_name=\${bin%.fa*}
        if [ ! -d \${bin_name}.tmp ]; then
            mkdir -p \${bin_name}.tmp
        fi

        bakta \\
            --db ${params.db} \\
            --threads ${task.cpus} \\
            --meta \\
            --output \${bin_name} \\
            --prefix \${bin_name} \\
            --tmp-dir \${bin_name}.tmp \\
            \${bin}

        if [ -d \${bin_name}.tmp ]; then
            rm -rf \${bin_name}.tmp
        fi
    done
    """

    stub:
    """
    for bin in ${bins}; do
        bin_name=\${bin%.fa*}
        mkdir \${bin_name}
        touch \${bin_name}/\${bin_name}.{embl,ffn,gbff,gff3,faa,fna}
        touch \${bin_name}/\${bin_name}.hypotheticals.{faa,tsv}
        touch \${bin_name}/\${bin_name}.{inference.tsv,png,svg,tsv,txt}
        touch \${bin_name}/\${bin_name}.log
    done
    """
}