#!/usr/bin/env python3
"""
Usage:
./consolidate_two_sets_of_bins.py \
    -1 <bin_folder_1> -2 <bin_folder_2> \
    -s1 <stats_file_1> -s2 <stats_file_2> \
    -o <output_folder> \
    -m min_completeness -x max_contamination
"""

import os
import argparse

from pathlib import Path


# Helper functions
def load_good_bins(stats_file: str, completeness: float, contamination: float) -> set:
    """ Load good bins from stats file

    Args:
        stats_file: str of Path to stats file

    Returns:
        good_bins: set of good bins
    """
    good_bins = set()
    for line in open(stats_file):
        if 'completeness' in line:
            continue
        cols = line.strip().split('\t')
        if float(cols[1]) >= completeness and float(cols[2]) <= contamination:
            good_bins.add(f'{cols[0]}.fa')


def load_contig_info(bin_file: Path, bins_set: set) -> dict:
    """ Load contig info from bin file

    Args:
        bin_file: Path to bin file
        bins_set: set of bins

    Returns:
        bins_info: dict of contig info
    """
    bins_info = {}
    for bin_file in bins_set:
        bins_info[bin_file] = {}
        contig_len = 0
        contig_name = ''
        for line in open(bin_file, 'rt'):
            if not line.startswith('>'):
                contig_len += len(line.strip())
            else:
                if contig_name:
                    bins_info[bin_file][contig_name] = contig_len
                    contig_len = 0
                contig_name = line.strip()[1:]

        bins_info[bin_file][contig_name] = contig_len

    return bins_info


def load_bins_stats(stats_file: Path) -> tuple:
    """ Load bins stats from stats file

    Args:
        stats_file: Path to stats file

    Returns:
        bins_stats, bins_summary: tuple of bins stats and summary
    """
    bins_stats, bins_summary = {}, {}
    for line in open(stats_file, 'rt'):
        if 'completeness' in line:
            bins_summary['header'] = line
            continue
        cols = line.strip().split('\t')
        bins_stats[f'{cols[0]}.fa'] = (float(cols[1]), float(cols[2]))
        bins_summary[f'{cols[0]}.fa'] = line

    return bins_stats, bins_summary


def existing_directory(path: str) -> Path:
    """ Check if path exists and is a directory

    If the path does not exist, an error will be raised.
    If the path is not a directory, an error will be raised.

    :param path: path string
    :return: Path object of absolute path
    """
    path = Path(path).absolute()
    if not path.exists():
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Not a directory: {path}")
    return path


def writable_directory(path: str) -> Path:
    """ Check if path exists and is a writable directory

    :param path: path string
    :return: Path object
    """
    path = Path(path).absolute()
    if path.exists():
        if not path.is_dir():
            raise argparse.ArgumentTypeError(f"Not a directory: {path}")
        if not os.access(path, os.W_OK):
            raise argparse.ArgumentTypeError(f"Directory not writable: {path}")
    return path


def parse_args() -> argparse.Namespace:
    """ Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Take two sets of bins and consolidate them')

    parser.add_argument('-1', '--binset-1', dest='bs1', metavar='BINSET1',
                        type=existing_directory, required=True,
                        help='Name of the first binsets directory')

    parser.add_argument('-2', '--binset-2', dest='bs2', metavar='BINSET2',
                        type=existing_directory, required=True,
                        help='Name of the second binsets directory')

    parser.add_argument('-s1', '--stats-1', dest='st1', metavar='STATS FILE 1',
                        type=str, required=True,
                        help='Stats file for the first binsets')

    parser.add_argument('-s2', '--stats-2', dest='st2', metavar='STATS FILE 2',
                        type=str, required=True,
                        help='Stats file for the second binsets')

    parser.add_argument('-o', '--output',
                        type=writable_directory, required=True,
                        help='Output directory')

    parser.add_argument('-m', '--min-completeness', type=float, default=50.,
                        help='Minimum completeness threshold [default: %(default)s]')
    parser.add_argument('-x', '--max-contamination', type=float, default=10.,
                        help='Maximum contamination threshold [default: %(default)s]')

    #parser.add_argument('-t', '--threads', type=int, default=1,
    #                    help='Number of threads [default: %(default)s]')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    output_dir = args.output
    comp = args.min_completeness
    cont = args.max_contamination

    good_bins_1 = load_good_bins(args.st1, comp, cont)
    good_bins_2 = load_good_bins(args.st2, comp, cont)

    # Load in the info about the contigs in each bin
    bins_1 = load_contig_info(args.bs1, good_bins_1)
    bins_2 = load_contig_info(args.bs2, good_bins_2)


    all_bin_pairs = {}
    for bin_1 in good_bins_1:
        all_bin_pairs[bin_1] = {}
        for bin_2 in good_bins_2:
            match_1_len = 0
            match_2_len = 0
            mismatch_1_len = 0
            mismatch_2_len = 0

            for contig in bins_1[bin_1]:
                if contig in bins_2[bin_2]:
                    match_1_len += bins_1[bin_1][contig]
                    match_2_len += bins_2[bin_2][contig]
                else:
                    mismatch_1_len += bins_1[bin_1][contig]
            for contig in bins_2[bin_2]:
                if contig in bins_1[bin_1]:
                    match_2_len += bins_2[bin_2][contig]
                else:
                    mismatch_2_len += bins_2[bin_2][contig]

            ratio_1 = 100 * match_1_len / (match_1_len + mismatch_1_len)
            ratio_2 = 100 * match_2_len / (match_2_len + mismatch_2_len)
            all_bin_pairs[bin_1][bin_2] = max(ratio_1, ratio_2)

    # Load in completeness and contamination scores of all the bins
    bins_1_stats, bins_1_summary = load_bins_stats(args.st1)
    bins_2_stats, bins_2_summary = load_bins_stats(args.st2)

    # Go through all good bins and choose the best ones
    output_dir.mkdir(parents=True, exist_ok=True)
    new_summary_file = bins_1_summary['header']
    bins_2_matches = {}
    bin_ct = 1
    for bin_1 in all_bin_pairs:
        score = bins_1_stats[bin_1][0] - bins_1_stats[bin_1][1] * 5
        found_better = False
        for bin_2 in all_bin_pairs[bin_1]:
            # check sufficient overlap (80% bin length)
            if all_bin_pairs[bin_1][bin_2] < 80:
                continue
            bins_2_matches[bin_2] = None
            # check if this bin is better than the current best
            new_score = bins_2_stats[bin_2][0] - bins_2_stats[bin_2][1] * 5
            if new_score > score:
                cmd = f'cp {args.bs2 / bin_2} {output_dir / f"bin.{bin_ct}.fa"}'
                new_summary_file += f'bin.{bin_ct}\t{'\t'.join(bins_2_summary[bin_2].split('\t')[1:])}'
                found_better = True
        if not found_better:
            new_summary_file += f'bin.{bin_ct}\t{'\t'.join(bins_1_summary[bin_1].split('\t')[1:])}'
            cmd = f'cp {args.bs1 / bin_1} {output_dir / f"bin.{bin_ct}.fa"}'
        os.system(cmd)
        bin_ct += 1

    print('Retrieve bins from second group that were not found in the first group')
    for bin_2 in bins_2_stats:
        if (bins_2_stats[bin_2][0] < comp) or (bins_2_stats[bin_2][1] > cont):
            continue
        if bin_2 in bins_2_matches:
            continue
        new_summary_file += f'bin.{bin_ct}\t{'\t'.join(bins_2_summary[bin_2].split('\t')[1:])}'
        cmd = f'cp {args.bs2 / bin_2} {output_dir / f"bin.{bin_ct}.fa"}'
        os.system(cmd)
        bin_ct += 1

    with open(output_dir / 'summary.txt', 'w') as f:
        f.write(new_summary_file)

    print(f'There were {bin_ct} bins cherry-picked from the original binsets.')