#!/usr/bin/env python3

import argparse
import os
import gzip

from pathlib import Path


def open_fasta(fasta_file):
    """ Open a fasta file

    Args:
        fasta_file (str): path to fasta file

    Returns:
        file: open file handle
    """

    if fasta_file.name.endswith('.gz'):
        return gzip.open(fasta_file, 'rt')
    else:
        return open(fasta_file, 'rt')


def read_contig_header(fasta_file):
    """ Read contig header

    Args:
        fasta_file (str): path to fasta file

    Yields:
        str: contig header
    """

    with open_fasta(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                yield line.strip().split()[0][1:]


def convert(paths, output_file: str, label: str = ''):
    """ Convert contigs to bin tsv for MAGScoT

    Args:
        paths (Path): path to binner output
        output_file (str): output file
        label (str, optional): label for binner. Defaults to ''.

    Returns:
        None
    """
    if label:
        label = f'\t{label}'

    fasta_exts = ('.fasta', '.fa', '.fna', '.bin', '_0.bin', '_1.bin', '.fa.gz')
    cur_binner = paths.name
    fasta_files = (file for file in paths.iterdir() if file.name.endswith(fasta_exts))

    with open(output_file, 'wt') as write_handler:
        for fasta_file in fasta_files:
            bin_id = fasta_file.stem
            for sequence_id in read_contig_header(fasta_file):
                write_handler.write(f"{cur_binner}_{bin_id}\t{sequence_id}{label}\n")


def existing_directory(path: str) -> Path:
    """ Check if path exists and is a directory

    If the path does not exist, an error will be raised.
    If the path is not a directory, an error will be raised.

    Args:
        path: path string

    Returns:
        Path object
    """
    path = Path(path).absolute()
    if not path.exists():
        raise argparse.ArgumentTypeError(f"Directory not found: {path}")
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"Not a directory: {path}")
    return path


def parse_args():
    """ Parse command line arguments

    Returns:
        argparse.Namespace: parsed arguments
    """

    parser = argparse.ArgumentParser(description='Convert contigs to bin tsv for MAGScoT')

    parser.add_argument('-p', '--paths', required=True, type=existing_directory,
                        help='Path(s) to binner output')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file')
    parser.add_argument('-l', '--binner-label', type=str, default='',
                        help='Label for binner')

    return parser.parse_args()


def main():
    args = parse_args()
    convert(args.paths, args.output, args.binner_label)


if __name__ == '__main__':
    main()