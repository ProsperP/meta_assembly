#!/usr/bin/env python3


import argparse
import os
import multiprocessing as mp

from pathlib import Path
from Bio import SeqIO


def run_multiprocess(task, iterable, cpus=1):
    """ Run a task in parallel

    Args:
        task: Function to run in parallel
        iterable: Iterable of arguments to pass to task
        cpus: Number of CPU cores to use

    Returns:
        List of results from task
    """
    cpus = min(mp.cpu_count(), cpus)
    with mp.Pool(processes=cpus) as pool:
        results = pool.map(func=task, iterable=iterable)
    return results


def create_dir(path: Path) -> None:
    """ Create a directory if it does not exist

    Args:
        path: Path object
    """
    if not path.exists():
        path.mkdir(parents=True)


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
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Rename metagenome-assembled genomes')

    parser.add_argument('-g', '--genome-dir', dest='input',
                        type=existing_directory, required=True,
                        help='Name of the genomes directory')

    parser.add_argument('-p', '--prefix', type=str, default='MAG',
                        help='MAGs ID prefix [default: %(default)s]')

    parser.add_argument('-o', '--output',
                        type=writable_directory, required=True,
                        help='Output directory')

    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads [default: %(default)s]')

    return parser.parse_args()


def main():
    args = parse_args()

    genome_dir = args.input
    bin_files = sorted([bin_file.name for bin_file in genome_dir.glob('*.fa')])

    outdir = args.output
    create_dir(outdir)

    new_mags_dir = outdir / 'MAGs_newid'
    if new_mags_dir.exists():
        print(f'Path `{new_mags_dir}` already existed. Would not override it. Exiting...')
        exit(1)
    else:
        create_dir(new_mags_dir)

    id_prefix = args.prefix
    mag_contig_sep = '-'
    id_convert_fh = (outdir / 'bins_mags_id.txt').open('wt')
    id_convert_fh.write('raw_id\tnew_id\n')
    bin_count = 1
    for bin_file in bin_files:
        old_bin = genome_dir / bin_file
        new_bin_id = f'{id_prefix}_{bin_count}'
        new_bin = new_mags_dir / f'{new_bin_id}.fa'
        bin_count += 1

        id_convert_fh.write(f'{bin_file.rstrip('.fa')}\t{new_bin_id}\n')

        with new_bin.open('wt') as outfh:
            contig_num = 1
            for record in SeqIO.parse(old_bin, 'fasta'):
                record.id = f'{new_bin_id}{mag_contig_sep}cg{contig_num}'
                contig_num += 1
                SeqIO.write(record, outfh, 'fasta')

    id_convert_fh.close()


if __name__ == '__main__':
    main()