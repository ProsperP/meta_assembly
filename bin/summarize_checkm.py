#!/usr/bin/env python3
""" Summarize checkM2 quality report

Contract (matched oringal metaWRAP usage):

1) summary_checkm.py -i checkm2_summary.tsv
   -> prints: bin   completeness    contamination   GC  lineage N50 size

2) summary_checkm.py -i checkm2_summary.tsv -l BINSET
   -> prints the same + trailing 'binner' column with BINSET for all rows

3) summary_checkm.py -i checkm2_summary.tsv -s binsM.stats
   -> ignores the -l arg (if set); uses binsM.stats as a source map
      key = bin name (col 1), value = binner label (col 8)
      If binsM.stats has a header and includes a 'binner' column, that
      column is used as the source of binner labels inseade col 8.
      Missing keys yield empty binner (no crash).
"""


import os
import sys
import argparse


# Define some helper functions
def fmtn(x, n: int = 5) -> str:
    """
    metaWRAP legacy 5-char -truncation.
    """
    return  '' if x is None else str(x)[:n]


def index_of(header: list, candidates: list) -> int:
    """
    Return the index of the first candidate in header (case-insensitive).
    If not found, return None.

    Args:
        header: list of str
        candidates: list of str

    Returns:
        int or None
    """
    lower_header = [h.strip().lower() for h in header]
    for word in candidates:
        cl = word.strip().lower()
        if cl in lower_header:
            return lower_header.index(cl)
    return None


def detect_indices(header: list) -> dict:
    """
    Detect the indices of required columns in checkM2 header.

    Args:
        header: list of str

    Returns:
        dict or None
    """
    # Build a dict of reuired column names from checkM2 header.
    # Robust to small naming variations.
    need = {
        'name': ['Name', 'Bin', 'bin', 'name'],
        'comp': ['Completeness', 'completeness', 'Completess (%)', 'completeness_%'],
        'cont': ['Contamination', 'contamination', 'Contamination (%)', 'contamination_%'],
        'gc':   ['GC_content', 'GC Content', 'GC', 'gc_content', 'gc'],
        'n50':  ['Contig_N50', 'N50', 'Contig n50'],
        'size': ['Genome_Size', 'Genome Size', 'Genome_size', 'Genome size', 'GenomeSize', 'Genome Size (bp)'],
    }
    idx = {}
    for key, candidates in need.items():
        i = index_of(header, candidates)
        if i is not None:
            idx[key] = i
        else:
            return None
    return idx


def load_source_map(path: str) -> dict[str, str]:
    """
    Return {bin_name: binner_label} from a checkM2 result TSV file.

    Primary behavior (orginal): use column 8 (index 7) as binner label.
    Graceful header handling: if header row exists and contains a 'binner'
    column, use that column's index instead.

    Args:
        path: a str path to a checkM2 result TSV file

    Returns:
        dict[str, str]
    """

    src = {}
    # Check if file exists
    if not os.path.isfile(path):
        return src

    with open(path, 'rt') as f:
        first_line = f.readline()
        if not first_line:
            return src

        columns = first_line.strip().split('\t')
        has_header = columns and (columns[0].strip().lower() == 'bin')
        binner_idx = 7
        if has_header:
            new_idx = index_of(columns, ['binner'])
            binner_idx = new_idx if new_idx else binner_idx
        else:
            src.update(parse_line(first_line, binner_idx))
        
        for line in f:
            src.update(parse_line(line, binner_idx))

    return src


def parse_line(line: str, binner_idx: int) -> dict:
    """
    Parse a line of checkM2 result TSV file.

    Args:
        line: a str line
        binner_idx: index of binner label column

    Returns:
        dict[str, str]
    """
    columns = line.strip().split('\t')
    if not columns:
        return {}

    name = columns[0]
    label = columns[binner_idx] if len(columns) > binner_idx else ''
    return {name: label}


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace
    """
    parser = argparse.ArgumentParser(description='Summarize checkM2 quality report')

    parser.add_argument('-i', '--report', dest='report', type=str, required=True,
                        metavar='quality_report.tsv', help='checkM2 result TSV file')
    parser.add_argument('-l', '--binner-label', dest='binner', type=str, help='binner label')
    parser.add_argument('-s', '--source-map', type=str, help='source map TSV file')
    parser.add_argument('-o', '--output', type=str, help='output file')

    return parser.parse_args()


def write_output(outfh, content: str) -> None:
    """
    Write content to output file or print to stdout.

    Args:
        outfh: file handler of output file
        content: str of content to be written

    Returns:
        None
    """
    if outfh:
        outfh.write(f'{content}\n')
    else:
        print(content)


def main():
    args = parse_args()

    base = ['bin', 'completeness', 'contamination', 'GC', 'lineage', 'N50', 'size']

    if args.output:
        outfh = open(args.output, 'wt')
    else:
        outfh = None

    if args.source_map:
        mode = 'map'
        binner_map = load_source_map(args.source_map)
    elif args.binner:
        mode = 'fixed'
        binner_fixed = sys.argv[2]
    else:
        # Default mode
        mode = 'plain'
        binner_fixed = None
        binner_map = None  

    with open(args.report, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        if header is None:
            sys.stderr.write(f'ERROR: empty quality report {inpath}\n')
            sys.exit(2)

        idx = detect_indices(header)
        if idx is None:
            sys.stderr.write(
                'ERROR: could not locate the required columns in'
                f' quality report header:\n {"\t".join(header)}\n'
            )
            sys.exit(2)

        if mode in ('map', 'fixed'):
            out_line = '\t'.join(base + ['binner'])
        else:
            out_line = '\t'.join(base)
        write_output(outfh, out_line)

        for row in reader:
            if not row: continue
            try:
                name = row[idx['name']]
            except IndexError:
                # malformed line; skip
                continue

            comp = row[idx['comp']] if idx['comp'] < len(row) else ''
            cont = row[idx['cont']] if idx['cont'] < len(row) else ''
            gc = row[idx['gc']] if idx['gc'] < len(row) else ''
            n50 = row[idx['n50']] if idx['n50'] < len(row) else ''
            size = row[idx['size']] if idx['size'] < len(row) else ''

            out = [name, fmtn(comp), fmtn(cont), fmtn(gc), 'NA', str(n50), str(size)]
            if mode == 'fixed':
                out.append(binner_fixed)
            elif mode == 'map':
                out.append(binner_map.get(name, ''))

            write_output(outfh, '\t'.join(out))

    if args.output:
        outfh.close()


if __name__ == '__main__':
    main()