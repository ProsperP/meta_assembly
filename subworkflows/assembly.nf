/*
 Assemble metagenomic genomes using metaSPAdes or MEGAHIT
 */

include {
    spades
    filterAssembly
    megahit
} from '../modules/local/assembly'


workflow ASSEMBLY {
    take:
    ch_clean_reads
    params

    main:

    ch_assembly_out = channel.empty()
    ch_assembly_log = channel.empty()

    if ( params.method == "spades" ) {
        spades(ch_clean_reads, params.spades)
        filterAssembly(spades.out.scaffolds, params.min_len)
        ch_assembly_out = ch_assembly_out.mix(
            spades.out.contigs,
            spades.out.scaffolds,
            spades.out.gene_clusters,
            spades.out.transcripts,
            spades.out.gfa
        ).ifEmpty([])
        ch_assembly_log = ch_assembly_log.mix(
            spades.out.warnings,
            spades.out.log
        )
        ch_final_assembly = filterAssembly.out.final_assembly
    } else if ( params.method == "megahit" ) {
        megahit(ch_clean_reads, params.min_len)
        ch_assembly_log = megahit.out.log
        ch_final_assembly = megahit.out.final_assembly
    }

    emit:
    final_assembly = ch_final_assembly
    assembly_out   = ch_assembly_out
    assembly_log   = ch_assembly_log
}