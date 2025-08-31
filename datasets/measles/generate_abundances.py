#!/usr/bin/env python3
"""Generate abundance table for Measles Penedos dataset.

Creates abundances.tsv with lines:
<filename>\t<abundance_int>

Same logic as Ebola generator (lognormal11 style): ceil( lognormal(mean, sigma) * scale ).
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
    ap = argparse.ArgumentParser(description='Generate abundances.tsv for Measles Penedos FASTA files (lognormal11 style).')
    ap.add_argument('--dir', type=Path, default=Path(__file__).parent, help='Directory containing FASTA(.gz) files (default: script directory)')
    ap.add_argument('--out', type=Path, default=None, help='Output TSV path (default: <dir>/abundances.tsv)')
    ap.add_argument('--mean', type=float, default=1.0, help='Mean parameter for lognormal (log-space)')
    ap.add_argument('--sigma', type=float, default=1.0, help='Sigma parameter for lognormal')
    ap.add_argument('--scale', type=float, default=10.0, help='Scale multiplier before integerization (default 10)')
    ap.add_argument('--seed', type=int, default=None, help='RNG seed for reproducibility')
    ap.add_argument('--rounding', choices=['ceil','round'], default='ceil', help='Integerization method')
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
        abundances.append(int(ceil(scaled)) if args.rounding == 'ceil' else int(round(scaled)))

    out_path = args.out if args.out is not None else directory / 'abundances.tsv'
    with open(out_path, 'w') as outf:
        for fname, abund in zip(fasta_files, abundances):
            outf.write(f'{fname}\t{abund}\n')
    total = sum(abundances)
    avg = total / len(abundances)
    print(f'Wrote {out_path} with {len(abundances)} entries. Total={total} Avg={avg:.2f} Seed={args.seed}')

if __name__ == '__main__':  # pragma: no cover
    main()
