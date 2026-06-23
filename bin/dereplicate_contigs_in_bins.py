#!/usr/bin/env python3
"""
Usage: dereplicate_contigs_in_bins.py bins.stats <bins_folder> <output_dir>
"""

import os
import sys


if __name__ == '__main__':
    # Load in bin completion and contamination scores
    bin_scores = {}
    for line in open(sys.argv[1], 'rt'):
        if 'completeness' in line:
            continue
        cols = line.strip().split('\t')
        score = float(cols[1]) - float(cols[2]) * 5 + 0.0000000001*float(cols[5])
        bin_scores[cols[0]] = score

    # Load in contigs in each bin
    contig_mapping = {}
    for bin_file in os.listdir(sys.argv[2]):
        bin_name = 
        for line in open(os.path.join(sys.argv[2], bin_file), 'rt'):
            if not line[0].startswith('>'):
                continue
            contig = line.strip()[1:]
            if contig not in contig_mapping:
                contig_mapping[contig] = bin_name
            else:
                if len(sys.argv) > 4:
                    if sys.argv[4] == 'remove':
                        contig_mapping[contig] = None
                    elif bin_scores[bin_name] > bin_scores[contig_mapping[contig]]:
                        contig_mapping[contig] = bin_name

    # Go over the bin files again and make a new dereplicated version of each bin file
    os.system(f'mkdir -p {sys.argv[3]}')
    for bin_file in os.listdir(sys.argv[2]):
        bin_name = 