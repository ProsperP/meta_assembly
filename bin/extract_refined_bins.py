#!/usr/bin/env python3
"""
Extract refined bins from a fasta file
"""

import argparse
import sys
import os
import gzip

from collections import defaultdict
from Bio import SeqIO


def open_fasta(fasta_file):
    """ Open a fasta file

    Args:
        fasta_file (str): path to fasta file

    Returns:
        file: open file handle
    """

    if fasta_file.endswith('.gz'):
        return gzip.open(fasta_file, 'rt')
    else:
        return open(fasta_file, 'rt')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("fasta_file", help="Input Fasta file.")
    parser.add_argument("cluster_file", help="Concoct output cluster file")
    parser.add_argument("--output_path", default=os.getcwd(), 
                        help="Directory where files will be printed")

    return parser.parse_args()


def main(args):
    all_seqs = {}
    for seq in SeqIO.parse(open_fasta(args.fasta_file), "fasta"):
        all_seqs[seq.id] = seq

    cluster_to_contigs = defaultdict(list)
    with open(args.cluster_file, 'rt') as fh:
        header = fh.readline().strip().split('\t')
        if header[0] != 'binnew' or header[1] != 'contig':
            sys.stderr.write("ERROR! Header line is not 'binnew, contig', please check your input file.\n")
            sys.exit(-1)

        for line in fh:
            line = line.strip().split('\t')
            cluster_to_contigs[line[0]].append(line[1])
    
    for bin_id, contig_ids in cluster_to_contigs.items():
        bin_id = bin_id.split('/')[-1]
        output_file = os.path.join(args.output_path, "{0}.fa".format(bin_id))
        seqs = [all_seqs[contig_id] for contig_id in contig_ids] 
        with open(output_file, 'wt') as ofh:
            SeqIO.write(seqs, ofh, 'fasta')


if __name__ == "__main__":
    args = parse_args()
    main(args)