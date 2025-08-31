import os, sys
from math import ceil
import gzip
from typing import List, Tuple

# Loader returns (genome_dir, genome_files, fixed_abundances_or_None)
# fixed_abundances_or_None is a list of floats (length = ngenomes) or None if distribution should be used later.

ALLOWED_FASTA_SUFFIXES = ('.fa', '.fna', '.fasta')
ALLOWED_FASTA_SUFFIXES_GZ = tuple(s + '.gz' for s in ALLOWED_FASTA_SUFFIXES) + ALLOWED_FASTA_SUFFIXES

def _list_fasta_files(directory: str) -> List[str]:
    return sorted([f for f in os.listdir(directory) if f.endswith(ALLOWED_FASTA_SUFFIXES_GZ)])

def _resolve_fname_with_gz(directory: str, fname: str) -> Tuple[str, str|None]:
    """Return existing filename (possibly switching to .gz) and optional info message."""
    base_path = os.path.join(directory, fname)
    if os.path.isfile(base_path):
        return fname, None
    # try gzip variant if not already gz
    if not fname.endswith('.gz') and os.path.isfile(base_path + '.gz'):
        return fname + '.gz', f'INFO: Using compressed file {fname}.gz instead of missing {fname}'
    # try removing .gz if provided
    if fname.endswith('.gz'):
        uncompressed = fname[:-3]
        if os.path.isfile(os.path.join(directory, uncompressed)):
            return uncompressed, f'INFO: Using uncompressed file {uncompressed} instead of missing {fname}'
    return fname, None  # will be caught by caller if truly missing

def load_dataset(name: str, ngenomes: int, user_abundances, allow_distribution: bool=True):
    """Load dataset metadata.

    name: one of ecoli, labmix, complex32, medium20, helicobacter-hepaticus, JGI
    ngenomes: number requested (must not exceed available)
    user_abundances: list[float] or None from -A
    allow_distribution: if False distribution will be ignored (when fixed abundances provided)

    Returns (genome_dir, genome_files, abundances_list_or_None, info_messages)
    info_messages: list of strings to print (INFO / WARNING messages)
    """
    msgs = []
    if name == 'ecoli':
        genome_dir = os.path.join('datasets','ecoli')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        genome_files = _list_fasta_files(genome_dir)
        if ngenomes > len(genome_files):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available ecoli genomes ({len(genome_files)})')
            sys.exit(1)
        genome_files = genome_files[:ngenomes]
        # If user provided abundances, we pass them back (validated in caller)
        return genome_dir, genome_files, user_abundances, msgs

    if name == 'labmix':
        genome_dir = os.path.join('datasets','labmix')
        abundance_file = os.path.join(genome_dir, 'abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for labmix dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available labmix genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for labmix; using abundances from abundances.tsv (truncated, no scaling).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'complex32':
        genome_dir = os.path.join('datasets','complex32')
        abundance_file = os.path.join(genome_dir, 'nanosim.abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for complex32 dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available complex32 genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        # Round abundances to nearest int as requested
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for complex32; using (rounded) abundances from nanosim.abundances.tsv (truncated).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'medium20':
        genome_dir = os.path.join('datasets','medium20')
        abundance_file = os.path.join(genome_dir, 'nanosim.abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for medium20 dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available medium20 genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for medium20; using (scaled x10, rounded) abundances from nanosim.abundances.tsv (truncated).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'JGI':
        genome_dir = os.path.join('datasets','JGI')  # directory name is capitalized
        abundance_file = os.path.join(genome_dir, 'nanosim.abundances.tsv')
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                if not os.path.isfile(os.path.join(genome_dir, actual_fname)):
                    msgs.append(f'WARNING: Listed genome missing in directory (even after .gz fallback): {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for JGI dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available JGI genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for JGI; using rounded abundances from nanosim.abundances.tsv (truncated).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'ebola':
        genome_dir = os.path.join('datasets','ebola')
        abundance_file = os.path.join(genome_dir, 'abundances.tsv')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}. Run generate_abundances.py first.')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                full_path = os.path.join(genome_dir, actual_fname)
                if not os.path.isfile(full_path):
                    msgs.append(f'WARNING: Listed genome missing in directory: {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for ebola dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available ebola genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for ebola; using abundances.tsv (already integer).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'measles':
        genome_dir = os.path.join('datasets','measles')
        abundance_file = os.path.join(genome_dir, 'abundances.tsv')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        if not os.path.isfile(abundance_file):
            print(f'ERROR: Missing abundance file {abundance_file}. Run generate_abundances.py first.')
            sys.exit(1)
        parsed = []
        with open(abundance_file) as af:
            for line in af:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                    continue
                fname, abund = parts
                try:
                    abund_val = float(abund)
                except ValueError:
                    msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                    continue
                actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                if info:
                    msgs.append(info)
                full_path = os.path.join(genome_dir, actual_fname)
                if not os.path.isfile(full_path):
                    msgs.append(f'WARNING: Listed genome missing in directory: {fname}')
                    continue
                parsed.append((actual_fname, abund_val))
        if len(parsed) == 0:
            print('ERROR: No valid genomes found for measles dataset.')
            sys.exit(1)
        if ngenomes > len(parsed):
            print(f'ERROR: --ngenomes ({ngenomes}) exceeds available measles genomes ({len(parsed)})')
            sys.exit(1)
        genome_files = [fname for fname, _ in parsed][:ngenomes]
        abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
        if user_abundances is not None:
            msgs.append('INFO: Ignoring --abundances for measles; using abundances.tsv (already integer).')
        return genome_dir, genome_files, abundances, msgs

    if name == 'sars_cov_2':
        genome_dir = os.path.join('datasets','sars_cov_2')
        if not os.path.isdir(genome_dir):
            print(f'ERROR: missing directory {genome_dir}')
            sys.exit(1)
        abundance_file = os.path.join(genome_dir, 'abundances.tsv')
        if os.path.isfile(abundance_file):
            parsed = []
            with open(abundance_file) as af:
                for line in af:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) != 2:
                        msgs.append(f'WARNING: Malformed line in abundances TSV ignored: {line}')
                        continue
                    fname, abund = parts
                    try:
                        abund_val = float(abund)
                    except ValueError:
                        msgs.append(f'WARNING: Non-numeric abundance for {fname}: {abund}; skipping')
                        continue
                    actual_fname, info = _resolve_fname_with_gz(genome_dir, fname)
                    if info:
                        msgs.append(info)
                    full_path = os.path.join(genome_dir, actual_fname)
                    if not os.path.isfile(full_path):
                        msgs.append(f'WARNING: Listed genome missing in directory: {fname}')
                        continue
                    parsed.append((actual_fname, abund_val))
            if len(parsed) == 0:
                print('ERROR: No valid genomes found for sars_cov_2 dataset (abundances.tsv present but empty/invalid).')
                sys.exit(1)
            if ngenomes > len(parsed):
                print(f'ERROR: --ngenomes ({ngenomes}) exceeds available sars_cov_2 genomes ({len(parsed)})')
                sys.exit(1)
            genome_files = [fname for fname, _ in parsed][:ngenomes]
            abundances = [int(round(abund)) for _, abund in parsed][:ngenomes]
            if user_abundances is not None:
                msgs.append('INFO: Ignoring --abundances for sars_cov_2; using abundances.tsv (already integer or rounded).')
            else:
                msgs.append('INFO: Using abundances from abundances.tsv for sars_cov_2.')
            return genome_dir, genome_files, abundances, msgs
        else:
            # Fallback: behave like ecoli (free distribution unless user supplies --abundances).
            genome_files = _list_fasta_files(genome_dir)
            if len(genome_files) == 0:
                print(f'ERROR: No FASTA(.gz) genomes found in {genome_dir}. Run download_sars_cov_2.py first.')
                sys.exit(1)
            if ngenomes > len(genome_files):
                print(f'ERROR: --ngenomes ({ngenomes}) exceeds available sars_cov_2 genomes ({len(genome_files)})')
                sys.exit(1)
            genome_files = genome_files[:ngenomes]
            if user_abundances is not None:
                msgs.append('INFO: Using explicit --abundances for sars_cov_2 (no abundances.tsv).')
            else:
                msgs.append('INFO: Using distribution sampling for sars_cov_2 (no abundances.tsv).')
            return genome_dir, genome_files, user_abundances, msgs

    print('ERROR: unknown dataset')
    sys.exit(1)
