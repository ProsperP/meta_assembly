#!/usr/bin/env python3
"""
Copyright (C) 2017, Weizhi Song, Torsten Thomas.
songwz03@gmail.com
t.thomas@unsw.edu.au

Binning_refiner is a free sofreware: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Binning_refiner is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

metaWRAP author notes:
I thank the original creator of this script! This is a great idea! To make
this script more usable as part of the metaWRAP binning pipeline, I removed
unnessary visual aspects of the original Bin_refiner script and made it
python2 compatible.

Check out the original program: https://github.com/songweizhi/Binning_refiner
And the publication: https://pubmed.ncbi.nlm.nih.gov/28186226/

modified by: Pu Jiajie <jiajie.pu@foxmail.com>
Author notes:
- Refactored code of the main() func into several helper functions.
- Added mulprocessing support.
"""


import os
import shutil
import gzip
import argparse

import multiprocessing as mp

from sys import stdout
from time import sleep
from pathlib import Path
from functools import partial

from Bio import SeqIO


# Define helper functions
def convert_to_mb(bp: int) -> float:
    """ Convert size in bytes to megabytes

    Args:
        size: Size in bytes
    
    Returns:
        Size in megabytes
    """
    return bp / (1024 * 1024)


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


def run_multiprocess_iter(task, iterable, cpus=1, chunksize=1):
    """Run a task in parallel and return an iterator (lazy evaluation).
    
    Args:
        task: Function to run in parallel
        iterable: Iterable of arguments to pass to task
        cpus: Number of CPU cores to use
        chunksize: Size of chunks to split the iterable into
    
    Yields:
        Results from task as they become available
    
    Example:
        >>> for result in run_multiprocess_iter(task, data, cpus=4):
        ...     print(result)
    """
    cpus = min(mp.cpu_count(), cpus)
    with mp.Pool(processes=cpus) as pool:
        yield from pool.imap(task, iterable, chunksize)


def check_bins_existence(folder: Path) -> tuple[bool, dict[Path, list[Path]]]:
    """ Check if a bin folder contains any bins

    Args:
        folder: Folder path to check for bin files
    
    Returns:
        Tuple containing:
            - bool: True if bins exist
            - dict: Mapping of folder to list of bin filenames
    """
    bins_files = [bin_file.name for bin_file in folder.glob('*.fa*')]
    return len(bins_files) > 0, {folder: bins_files}


def has_same_extension(files: list) -> bool:
    """ Check if all files have the same extension

    Args:
        files: a list of bin file names

    Returns:
        bool: True if all files have the same extension
    """
    file_ext_set = set(Path(fname).suffix for fname in files)
    return len(file_ext_set) == 1


def validate_and_collect_bins(res_list: list[Path]) -> dict[str, list[str]]:
    """ Validate and collect bin files

    Args:
        res_list: list of results from check_bins_existence

    Returns:
        dict[str, list[str]]: Mapping of folder to list of bin filenames
    """
    folder_bins_dict = {}
    for has_bin, res_dict in res_list:
        folder = list(res_dict)[0]
        if not has_bin:
            print(f'No bins found in {folder}, please double check!')
            exit()
        if not has_same_extension(res_dict[folder]):
            print(f'Bins in {folder} have different file extensions, please use '
                  'same extension (fa, fas or fasta) for all bins in the same folder.')
            exit()
        folder_bins_dict.update(res_dict)
    return folder_bins_dict


def create_dir(path: Path) -> None:
    """ Create a directory if it does not exist

    Args:
        path: Path object
    """
    if not path.exists():
        path.mkdir(parents=True)


def get_open_function(file_name: str):
    """ Get the appropriate open function based on file extension

    Args:
        file_name: File name or path

    Returns:
        gzip.open if file ends with .gz, otherwise open
    """
    if file_name.endswith('.gz'):
        return gzip.open
    return open


def prepend_bin_info_to_contigs(
        bin_file: str,
        parent_dir: Path,
        open_func,
        temp_dir: Path,
        sep: str = '__') -> None:
    """ Create folder to hold bins with renamed contig name

    Args:
        bin_file: Bin file name
        parent_dir: Parent directory of bin files
        open_func: Open function
        temp_dir: Temporary directory to store renamed bins
        sep: Separator to use in new contig name

    Returns:
        None
    """
    file_path = Path(parent_dir, bin_file)
    fname = bin_file.rstrip(''.join(file_path.suffixes))
    infh = open_func(file_path, 'rt')
    bin_content = SeqIO.parse(infh, 'fasta')
    new_file = Path(temp_dir, f'{parent_dir.name}_{fname}.fasta').open('wt')
    for contig in bin_content:
        new_id = f'{parent_dir.name}{sep}{fname}{sep}{contig.id}'
        contig.id = new_id
        contig.description = ''
        SeqIO.write(contig, new_file, 'fasta')
    infh.close()
    new_file.close()


def write_consensus_contigs(file_name: str, contig_bin_dict: dict, binner_num: int):
    """ Write consensus contigs to file
    
    Args:
        file_name: File name to write consensus contigs
        contig_bin_dict: Dictionary of contig to bins
        binner_num: Number of bins

    Returns:
        None
    """
    with open(file_name, 'wt') as contig_assignment:
        for record_id in contig_bin_dict:
            # if all bins have this contig
            if len(contig_bin_dict[record_id]) == binner_num:
                contig_id, contig_len = record_id
                contig_assignment.write(
                    '{}\t{}\t{}\n'
                    .format('\t'.join(contig_bin_dict[record_id]),
                            contig_id, contig_len)
                )


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


def writable_directory(path: str) -> Path:
    """ Check if path exists and is a writable directory

    Args:
        path: path string

    Returns:
        Path object
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
    parser = argparse.ArgumentParser(description='Refine binning results')

    parser.add_argument('-1', '--binset-1', dest='bs1', metavar='BINSET1',
                        type=existing_directory, required=True,
                        help='Name of the first binsets directory')

    parser.add_argument('-2', '--binset-2', dest='bs2', metavar='BINSET2',
                        type=existing_directory, required=True,
                        help='Name of the second binsets directory')

    parser.add_argument('-3', '--binset-3', dest='bs3', metavar='BINSET3',
                        type=existing_directory, required=False,
                        help='Name of the third binsets directory')

    parser.add_argument('-ms', '--min-size', type=int, default=524288,
                        help='Minimum size for refined bins [default: %(default)s]')

    parser.add_argument('-o', '--output',
                        type=writable_directory, required=True,
                        help='Output directory')

    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads [default: %(default)s]')

    return parser.parse_args()


def main():
    args = parse_args()

    bin_size_cutoff = args.min_size
    bin_size_cutoff_mb = convert_to_mb(bin_size_cutoff)

    outdir = args.output

    bin_folders = [args.bs1, args.bs2]
    if args.bs3:
        bin_folders.append(args.bs3)

    # Validate and collect input bins
    res_list = run_multiprocess(check_bins_existence, bin_folders,
                                cpus=len(bin_folders))
    folder_bins_dict = validate_and_collect_bins(res_list)

    # Create output directory after input validation
    create_dir(outdir)

    # Create folder to hold bins with renamed contig name
    sep = '__'
    tempdirs = []
    for each_folder in bin_folders:
        # TODO: use logging instead of print
        print(f'Add folder/bin name to config name for bins in {each_folder}')
        temp_dir = outdir.joinpath(f'{each_folder.name}_temp')
        tempdirs.append(temp_dir)
        create_dir(temp_dir)

        # Add binning program and bin id to bin's contig name
        bin_files = folder_bins_dict[each_folder]
        open_func = get_open_function(bin_files[0])
        # Prepare for multiprocessing
        partial_func = partial(prepend_bin_info_to_contigs,
                               parent_dir=each_folder,
                               open_func=open_func,
                               temp_dir=temp_dir, sep=sep)
        run_multiprocess(partial_func, bin_files, cpus=args.threads)

    # Combine all modified bins into one file
    temp_fas = [str(temp_dir.joinpath('*.fasta')) for temp_dir in tempdirs]
    # Save for later use
    input_contigs_file = outdir.joinpath(f'combined_{bin_folders[0].name}_bins.fasta')
    os.system(f'cat {temp_fas[0]} > {input_contigs_file}')

    combined_bins_file = outdir.joinpath('combined_all_bins.fasta')
    os.system(f'cat {' '.join(temp_fas)} > {combined_bins_file}')
    # Clean up
    os.system(f'rm -r {' '.join(str(temp_dir) for temp_dir in tempdirs)}')


    # TODO: write into a function
    combined_bins = SeqIO.parse(combined_bins_file, 'fasta')
    contig_bin_dict = {}
    for contig in combined_bins:
        id_split = contig.id.split(sep)
        dir_name = id_split[0]
        bin_name = id_split[1]
        # assign record id as (contig_id, contig_len)
        record_id = (id_split[2], len(contig.seq))

        if record_id not in contig_bin_dict:
            contig_bin_dict[record_id] = [f'{dir_name}{sep}{bin_name}']
        else:
            contig_bin_dict[record_id].append(f'{dir_name}{sep}{bin_name}')

    contig_assign_file = outdir.joinpath('contig_assignment.txt')
    # Only keep contigs existed in all bin sets!
    write_consensus_contigs(contig_assign_file, contig_bin_dict, len(bin_folders))

    contig_assign_file_sort = outdir.joinpath('contig_assignment_sort.txt')
    contig_assign_file_sort_one_line = outdir.joinpath('contig_assignment_sort_one_line.txt')
    os.system(f'sort {contig_assign_file} > {contig_assign_file_sort}')

    # TODO: write into a function
    contig_assign_sort = open(contig_assign_file_sort, 'rt')
    contig_assign_sort_one_line = open(contig_assign_file_sort_one_line, 'wt')
    current_match = None
    current_match_contigs = []
    current_len_total = 0
    qualified_bins_num = 0  # Total number of qualified bins
    for line in contig_assign_sort:
        cols = line.strip().split('\t')
        current_contig = cols[2]
        current_len = int(cols[3])
        matched_bins = '\t'.join(cols[:2])
        if current_match is None:
            current_match = matched_bins
            current_match_contigs.append(current_contig)
            current_len_total += current_len
        elif current_match == matched_bins:
            current_match_contigs.append(current_contig)
            current_len_total += current_len
        elif current_match != matched_bins:
            if current_len_total >= bin_size_cutoff:
                qualified_bins_num += 1
                contig_assign_sort_one_line.write(
                    'Refined_{}\t{}\t{}\t{}\n'
                    .format(qualified_bins_num, current_match, current_len_total,
                            '\t'.join(current_match_contigs))
                )

            current_match = matched_bins
            current_match_contigs = [current_contig]
            current_len_total = current_len

    if current_len_total >= bin_size_cutoff:
        qualified_bins_num += 1
        contig_assign_sort_one_line.write(
            'Refined_{}\t{}\t{}\t{}\n'
            .format(qualified_bins_num, current_match, current_len_total,
                    '\t'.join(current_match_contigs))
        )
    contig_assign_sort_one_line.close()

    print(f'The number of refined bins: {qualified_bins_num}')


    # TODO: export refined bins and prepare input for GoogleVis
    print('Exporting refined bins...')
    bin_src_len_file = outdir.joinpath('Refined_bins_sources_and_length.txt')
    bin_contig_file = outdir.joinpath('Refined_bins_contigs.txt')
    googlevis_file = outdir.joinpath(f'GoogleVis_Sankey_{bin_size_cutoff_mb:.2f}Mbp.csv')
    refined_outdir = outdir.joinpath('Refined_bins')
    create_dir(refined_outdir)

    refined_bins = open(contig_assign_file_sort_one_line, 'rt')
    bin_src_len_fh = open(bin_src_len_file, 'wt')
    bin_contig_fh = open(bin_contig_file, 'wt')
    googlevis_fh = open(googlevis_file, 'wt')

    googlevis_fh.write('C1,C2,Length (Mbp)\n')
    raw_bin_dir_num = len(bin_folders)
    for refined_bin in refined_bins:
        cols = refined_bin.strip().split('\t')
        bin_name = cols[0]
        bin_src = cols[1:raw_bin_dir_num+1]
        bin_len = int(cols[raw_bin_dir_num+1])
        bin_contig = cols[raw_bin_dir_num+2:]
        bin_src_len_fh.write(f'{bin_name}\t{bin_len}bp\t{'\t'.join(bin_src)}\n')
        bin_contig_fh.write(f'{bin_name}\n{"\t".join(bin_contig)}\n')
        bin_len_mbp = convert_to_mb(bin_len)

        for i in range(len(bin_src)-1):
            googlevis_fh.write(f'{bin_src[i]},{bin_src[i+1]},{bin_len_mbp:.2f}\n')

        print(f'\rExtracting refined bins: {bin_name}.fasta', end='')
        refined_bin_file = refined_outdir.joinpath(f'{bin_name}.fasta')
        refined_bin_fh = open(refined_bin_file, 'wt')
        input_contigs = SeqIO.parse(input_contigs_file, 'fasta')
        for input_contig in input_contigs:
            input_contig_id = input_contig.id.split(sep)[-1]
            if input_contig_id in bin_contig:
                input_contig.id = input_contig_id
                input_contig.description = ''
                SeqIO.write(input_contig, refined_bin_fh, 'fasta')
        refined_bin_fh.close()

    googlevis_fh.close()
    bin_src_len_fh.close()
    bin_contig_fh.close()
    
    # Remove temp files
    print('\nCleaning up...')
    os.system(f'rm {contig_assign_file} {contig_assign_file_sort}')
    os.system(f'rm {contig_assign_file_sort_one_line}')
    os.system(f'rm {combined_bins_file} {input_contigs_file}')
    print('\nAll done!')


if __name__ == '__main__':
    main()