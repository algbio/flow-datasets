#!/usr/bin/env python3
"""Generate abundance table for Ebola Makona dataset.

Creates abundances.tsv in the current directory with lines:
<filename>\t<abundance_int>

Abundances are sampled using the same logic as the "lognormal11" distribution in
construct-flow-graph.py:
    ceil( lognormal(mean=1, sigma=1) * 10 )

Options let you override mean, sigma, scale multiplier, RNG seed, and whether to
use ceil vs round.

Basic usage:
    # From repository root (creates datasets/ebola/abundances.tsv)
    python datasets/ebola/generate_abundances.py --seed 42

Specify a different output path and parameters:
    python datasets/ebola/generate_abundances.py \
        --dir datasets/ebola \
        --out datasets/ebola/abundances.custom.tsv \
        --mean 1.0 --sigma 1.0 --scale 10 --seed 1234

Result format (TSV):
    <filename>\t<integer_abundance>
"""
from __future__ import annotations
from pathlib import Path
from math import ceil
import argparse
from numpy.random import default_rng

ALLOWED_SUFFIXES = ('.fa', '.fna', '.fasta')
ALLOWED_SUFFIXES_GZ = tuple(s + '.gz' for s in ALLOWED_SUFFIXES) + ALLOWED_SUFFIXES

def list_fasta(directory: Path):
    return sorted([p.name for p in directory.iterdir() if p.is_file() and p.name.endswith(ALLOWED_SUFFIXES_GZ)])

def main():
    ap = argparse.ArgumentParser(description="Generate abundances.tsv for Ebola Makona FASTA files (lognormal11 style).")
    ap.add_argument('--dir', type=Path, default=Path(__file__).parent, help='Directory containing FASTA(.gz) files (default: script directory)')
    ap.add_argument('--out', type=Path, default=None, help='Output TSV path (default: <dir>/abundances.tsv)')
    ap.add_argument('--mean', type=float, default=1.0, help='Mean parameter for lognormal (in log-space)')
    ap.add_argument('--sigma', type=float, default=1.0, help='Sigma parameter for lognormal')
    ap.add_argument('--scale', type=float, default=10.0, help='Scale multiplier before integerization (default 10)')
    ap.add_argument('--seed', type=int, default=None, help='RNG seed for reproducibility')
    ap.add_argument('--rounding', choices=['ceil','round'], default='ceil', help='Integerization method (ceil mimics existing lognormal11 logic)')
    args = ap.parse_args()

    rng = default_rng(args.seed)
    directory = args.dir.resolve()
    if not directory.is_dir():
        raise SystemExit(f'ERROR: directory not found: {directory}')

    fasta_files = list_fasta(directory)
    if not fasta_files:
        raise SystemExit(f'ERROR: no FASTA(.gz) files found in {directory}')

    samples = rng.lognormal(mean=args.mean, sigma=args.sigma, size=len(fasta_files))
    abundances = []
    for x in samples:
        scaled = x * args.scale
        if args.rounding == 'ceil':
            abundances.append(int(ceil(scaled)))
        else:
            abundances.append(int(round(scaled)))

    out_path = args.out if args.out is not None else directory / 'abundances.tsv'
    with open(out_path, 'w') as outf:
        for fname, abund in zip(fasta_files, abundances):
            outf.write(f'{fname}\t{abund}\n')

    # Basic summary to stderr/stdout
    total = sum(abundances)
    avg = total / len(abundances)
    print(f'Wrote {out_path} with {len(abundances)} entries. Total={total} Avg={avg:.2f} Seed={args.seed}')

if __name__ == '__main__':
    main()
