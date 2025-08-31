#!/usr/bin/env python3
"""Download Measles virus complete genomes.

Primary mode (paper-specific): attempt to discover accessions linked to
Penedos et al. 2015 (PMID 26569100) via PubMed -> nuccore elink.

If that yields none (common because not all measles submissions are cross-linked),
fallback strategies are tried in order:
    1. ESearch query (default configurable) against nuccore for complete measles genomes.
    2. Static curated accession list (ship a modest deterministic subset).

Each retained sequence must have a FASTA header containing both "complete genome" and
"Measles" (case-insensitive) and length between 15,000 and 17,000 nt (basic sanity
for MeV ~15.8 kb). The extra "Measles" gate reduces accidental inclusion of
unrelated sequences that match a generic complete genome query.
One compressed FASTA per accession: <accession>.fasta.gz
"""
from __future__ import annotations
from pathlib import Path
from time import sleep
from Bio import Entrez
from xml.etree import ElementTree
import ssl
import hashlib
import certifi
import urllib.request
import urllib.error
import argparse
import sys
import gzip

HERE = Path(__file__).parent
PMID = '26569100'  # Penedos et al. 2015 (paper-specific attempt)

# Parameters
DEFAULT_EMAIL = 'flow-datasets@example.com'  # fallback if user does not configure
Entrez.tool = 'create-flow-graphs-downloader'

def configure_entrez(email: str):
    Entrez.email = email
    # Create an SSL context using certifi to avoid local certificate store issues.
    _SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
    def _entrez_urlopen(*args, **kwargs):  # type: ignore[no-untyped-def]
        if 'context' not in kwargs:
            kwargs['context'] = _SSL_CONTEXT
        return urllib.request.urlopen(*args, **kwargs)
    Entrez.urlopen = _entrez_urlopen  # type: ignore[attr-defined]


def fetch_fasta(acc: str, retries: int = 3, delay: float = 0.5) -> str:
    for attempt in range(1, retries + 1):
        try:
            with Entrez.efetch(db='nuccore', id=acc, rettype='fasta', retmode='text') as handle:
                return handle.read()
        except urllib.error.URLError as e:
            if attempt == retries:
                raise
            print(f'Retry {attempt}/{retries} after URLError for {acc}: {e.reason}', file=sys.stderr)
            sleep(delay * attempt)
        except Exception as e:  # noqa: BLE001
            if attempt == retries:
                raise
            print(f'Retry {attempt}/{retries} after error for {acc}: {e}', file=sys.stderr)
            sleep(delay * attempt)
    return ''


def fetch_all_accessions_from_pubmed(pmid: str) -> list[str]:
    with Entrez.elink(dbfrom='pubmed', id=pmid, linkname='pubmed_nuccore') as h:
        link_root = ElementTree.parse(h).getroot()
    nuccore_ids = [e.text for e in link_root.findall('.//LinkSetDb/Link/Id')]
    if not nuccore_ids:
        return []
    accs: list[str] = []
    for i in range(0, len(nuccore_ids), 200):
        batch = nuccore_ids[i:i+200]
        with Entrez.esummary(db='nuccore', id=','.join(batch)) as h:
            summ_root = ElementTree.parse(h).getroot()
        for doc in summ_root.findall('.//DocSum'):
            caption = None
            version = None
            for item in doc.findall('Item'):
                if item.get('Name') == 'Caption':
                    caption = item.text
                elif item.get('Name') == 'Version':
                    version = item.text
            if caption and version:
                accs.append(f"{caption}.{version}")
    return sorted(set(accs))

# Static curated fallback (first page + a few classic strains) to guarantee reproducibility.
# You can extend this list if you need more diversity.
STATIC_FALLBACK: list[str] = [
    "K01711.1", "AF266288.1", "AF266291.1", "AF266290.1", "PX121934.1", "PX121931.1", "PV931737.1", "PV931736.1", "PV870366.1", "PV870365.1",
    "PV870364.1", "PV870363.1", "PV870362.1", "PV870361.1", "PV870360.1", "PV870359.1", "PV870358.1", "PV870357.1", "PV870356.1", "PV870355.1",
    "PV870354.1", "PV870353.1", "PV870352.1", "JN635407.1", "JN635408.1", "JN635409.1", "KT732231.1", "KT732232.1", "KT732233.1", "MH638233.1",
]

def esearch_accessions(query: str, retmax: int = 100) -> list[str]:
    """Return accession.version list for an arbitrary nuccore esearch query.

    We map IDs -> summaries to extract Caption+Version (accession.version).
    """
    with Entrez.esearch(db='nuccore', term=query, retmax=retmax) as h:
        rec = Entrez.read(h)
    id_list = rec.get('IdList', [])
    accs: list[str] = []
    if not id_list:
        return accs
    for i in range(0, len(id_list), 200):
        batch = id_list[i:i+200]
        with Entrez.esummary(db='nuccore', id=','.join(batch)) as h:
            summ_root = ElementTree.parse(h).getroot()
        for doc in summ_root.findall('.//DocSum'):
            caption = None
            version = None
            length = None
            for item in doc.findall('Item'):
                name = item.get('Name')
                if name == 'Caption':
                    caption = item.text
                elif name == 'Version':
                    version = item.text
                elif name == 'Length':
                    length = item.text
            if caption and version:
                if length is not None:
                    try:
                        L = int(length)
                        if not (15000 <= L <= 17000):  # crude measles genome length gate
                            continue
                    except ValueError:
                        pass
                accs.append(f'{caption}.{version}')
    return sorted(set(accs))


def main():
    ap = argparse.ArgumentParser(description='Download Measles virus complete genomes (paper-linked or generic).')
    # Email removed from required args; use default (NCBI recommends providing one – adjust in code if desired)
    ap.add_argument('--sleep', type=float, default=0.34, help='Delay between requests (~3/sec).')
    ap.add_argument('--max', type=int, default=None, help='Limit number of accessions (after selection).')
    ap.add_argument('--query', default='"Measles morbillivirus"[Organism] AND "complete genome"',
                    help='ESearch fallback query used if PubMed elink yields none.')
    ap.add_argument('--retmax', type=int, default=120, help='Retmax for esearch fallback.')
    ap.add_argument('--no-elink', action='store_true', help='Skip PubMed elink step and go straight to esearch.')
    ap.add_argument('--validate-static', action='store_true',
                    help='Validate STATIC_FALLBACK accessions (existence, header phrase, length) and output a report; does not perform normal discovery logic.')
    ap.add_argument('--report', default='static_validation.tsv',
                    help='Report filename for --validate-static (TSV).')
    ap.add_argument('--emit-python-list', action='store_true',
                    help='With --validate-static, also emit a cleaned Python list of valid accessions to stdout at end.')
    args = ap.parse_args()

    # Configure Entrez with built-in default email (edit DEFAULT_EMAIL in source to customize)
    configure_entrez(DEFAULT_EMAIL)

    # Optional validation-only mode
    if args.validate_static:
        print(f'Validating {len(STATIC_FALLBACK)} static fallback accessions...')
        report_path = HERE / args.report
        valid: list[str] = []
        with open(report_path, 'w') as rep:
            rep.write('accession\tstatus\tlength\theader_has_complete\theader_has_measles\treason\n')
            for i, acc in enumerate(STATIC_FALLBACK, 1):
                try:
                    seq_txt = fetch_fasta(acc)
                    if not seq_txt:
                        rep.write(f'{acc}\tfail\tNA\tNA\tempty response\n')
                        print(f'{acc}: FAIL empty')
                        continue
                    header = seq_txt.splitlines()[0]
                    header_lc = header.lower()
                    has_complete = 'complete genome' in header_lc
                    has_measles = 'measles' in header_lc
                    seq_lines = [l.strip() for l in seq_txt.splitlines()[1:] if l and not l.startswith('>')]
                    seq_len = sum(len(l) for l in seq_lines)
                    if not has_complete:
                        rep.write(f'{acc}\tfail\t{seq_len}\t{int(has_complete)}\t{int(has_measles)}\tmissing "complete genome" phrase\n')
                        print(f'{acc}: FAIL missing "complete genome" (len {seq_len})')
                        continue
                    if not has_measles:
                        rep.write(f'{acc}\tfail\t{seq_len}\t{int(has_complete)}\t{int(has_measles)}\tmissing "Measles" token\n')
                        print(f'{acc}: FAIL missing "Measles" token (len {seq_len})')
                        continue
                    if not (15000 <= seq_len <= 16500):
                        rep.write(f'{acc}\tfail\t{seq_len}\t{int(has_complete)}\t{int(has_measles)}\tlength out of range\n')
                        print(f'{acc}: FAIL length {seq_len}')
                        continue
                    # Passed all gates
                    rep.write(f'{acc}\tok\t{seq_len}\t{int(has_complete)}\t{int(has_measles)}\t\n')
                    valid.append(acc)
                    print(f'{acc}: OK (len {seq_len})')
                except Exception as e:  # noqa: BLE001
                    rep.write(f'{acc}\tfail\tNA\tNA\tNA\t{e.__class__.__name__}:{e}\n')
                    print(f'{acc}: FAIL exception {e}')
                sleep(args.sleep)
        print(f'Validation complete. {len(valid)}/{len(STATIC_FALLBACK)} valid. Report: {report_path}')
        if args.emit_python_list and valid:
            cleaned = ','.join(f'"{v}"' for v in valid)
            print('\n# Cleaned STATIC_FALLBACK (valid only)\n[')
            # chunk output for readability
            for j in range(0, len(valid), 10):
                chunk = valid[j:j+10]
                print('    ' + ', '.join(f'"{c}"' for c in chunk) + ',')
            print(']')
        return 0

    accessions: list[str] = []
    if not args.no_elink:
        print('Attempting PubMed -> nuccore elink (PMID 26569100)...')
        try:
            accessions = fetch_all_accessions_from_pubmed(PMID)
        except Exception as e:  # noqa: BLE001
            print(f'Warning: elink retrieval failed: {e}')
    if not accessions:
        print('No (or zero) elink accessions; attempting ESearch fallback...')
        try:
            accessions = esearch_accessions(args.query, retmax=args.retmax)
        except Exception as e:  # noqa: BLE001
            print(f'ESearch fallback failed: {e}')
    if not accessions:
        print(f'Both elink and esearch produced no accessions; using static curated fallback ({len(STATIC_FALLBACK)}).')
        accessions = STATIC_FALLBACK
    if not accessions:
        print('ERROR: Still no accessions available. Aborting.')
        return 1

    # Stable ordering
    accessions = sorted(set(accessions))
    if args.max is not None:
        accessions = accessions[:args.max]
    print(f'Proceeding with {len(accessions)} candidate accessions.')

    kept = 0
    duplicate_skipped = 0
    seen_hashes: dict[str, str] = {}
    for i, acc in enumerate(accessions, 1):
        try:
            seq_txt = fetch_fasta(acc)
            if not seq_txt:
                print(f'Empty response for {acc}')
                continue
            header = seq_txt.splitlines()[0] if seq_txt else ''
            header_lc = header.lower()
            if 'complete genome' not in header_lc:
                print(f'Skipping {acc}: header lacks "complete genome"')
                continue
            if 'measles' not in header_lc:
                print(f'Skipping {acc}: header lacks "Measles"')
                continue
            # Optional length gate (cheap) from sequence itself
            seq_lines = [l.strip() for l in seq_txt.splitlines()[1:] if l and not l.startswith('>')]
            seq_len = sum(len(l) for l in seq_lines)
            if not (15000 <= seq_len <= 16500):
                print(f'Skipping {acc}: length {seq_len} outside 15k-16.5k window')
                continue
            # Deduplicate identical genome sequences across different accessions (default behavior)
            seq_only = ''.join(seq_lines).upper()
            h = hashlib.sha256(seq_only.encode('utf-8')).hexdigest()
            prev = seen_hashes.get(h)
            if prev is not None:
                print(f'Skipping {acc}: identical sequence to {prev}')
                duplicate_skipped += 1
                continue
            out_path = HERE / f'{acc}.fasta.gz'
            if out_path.exists():
                print(f'[exists] {out_path.name}')
            else:
                with gzip.open(out_path, 'wt') as gz:
                    gz.write(seq_txt if seq_txt.startswith('>') else f'>{acc}\n{seq_txt}')
                kept += 1
                seen_hashes[h] = acc
        except Exception as e:  # noqa: BLE001
            print(f'Failed {acc}: {e}')
        if i % 5 == 0 or i == len(accessions):
            print(f'Processed {i}/{len(accessions)} (new kept {kept})')
        sleep(args.sleep)

    print(f'Finished. Downloaded {kept} unique complete genome sequences (duplicates skipped: {duplicate_skipped}).')
    return 0

if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
