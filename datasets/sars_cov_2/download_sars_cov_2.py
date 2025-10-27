#!/usr/bin/env python3
"""Download up to 100 SARS-CoV-2 complete genome sequences from NCBI nuccore.

Strategy:
- Use an Entrez esearch query restricted to SARS-CoV-2 organism name and "complete genome" in the title/definition, with a broad publication date window.
- Retrieve accession.version via esummary (preferring AccessionVersion; fallback Caption+Version) for the result UIDs.
- Filter sequences to expected genome length (29k-31k) and headers containing 'complete genome'.
- Additional quality filter: retain only sequences whose fraction of 'N' (ambiguous) bases is < 10% (case-insensitive) across the concatenated sequence lines.
- Deterministically select the first 100 unique accession versions after lexicographic sort.
- Download each as its own compressed FASTA <accession>.fasta.gz.

Reproducibility: the chosen accession list is written to both selected_accessions.txt and
static_accessions.txt (the latter can be committed and reused to skip dynamic querying).
"""
from __future__ import annotations
from pathlib import Path
from time import sleep
from Bio import Entrez
from xml.etree import ElementTree
import ssl
import certifi
import urllib.request
import urllib.error
import gzip
import sys

HERE = Path(__file__).parent
MAX_ACCESSIONS = 100
MAX_N_FRACTION = 0.10  # keep sequence only if (N count / sequence length) < this threshold

Entrez.email = 'flow-datasets@example.com'  # set your email if desired
Entrez.tool = 'create-flow-graphs-downloader'

_SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())

def _entrez_urlopen(*args, **kwargs):  # type: ignore[no-untyped-def]
    if 'context' not in kwargs:
        kwargs['context'] = _SSL_CONTEXT
    return urllib.request.urlopen(*args, **kwargs)

Entrez.urlopen = _entrez_urlopen  # type: ignore[attr-defined]

QUERY = '"Severe acute respiratory syndrome coronavirus 2"[Organism] AND complete genome[Title] AND ("2020/01/01"[PDAT] : "3000"[PDAT])'


def esearch_uids(query: str) -> list[str]:
    with Entrez.esearch(db='nuccore', term=query, retmax=5000) as h:
        root = ElementTree.parse(h).getroot()
    return [e.text for e in root.findall('.//IdList/Id')]


def esummary_accessions(uids: list[str]) -> list[tuple[str, int, str]]:
    accs: list[tuple[str, int, str]] = []  # (acc.version, length, title)
    for i in range(0, len(uids), 200):
        batch = uids[i:i+200]
        with Entrez.esummary(db='nuccore', id=','.join(batch)) as h:
            sroot = ElementTree.parse(h).getroot()
        for doc in sroot.findall('.//DocSum'):
            data = {}
            for item in doc.findall('Item'):
                data[item.get('Name')] = item.text
            acc_ver = data.get('AccessionVersion')
            if not acc_ver and data.get('Caption') and data.get('Version'):
                acc_ver = f"{data.get('Caption')}.{data.get('Version')}"
            length = data.get('Length')
            title = data.get('Title')
            if acc_ver and length and title:
                try:
                    ln = int(length)
                except ValueError:
                    continue
                accs.append((acc_ver, ln, title))
    return accs


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


def main():  # noqa: D401
    print('Searching for SARS-CoV-2 complete genomes...')
    try:
        uids = esearch_uids(QUERY)
    except Exception as e:  # noqa: BLE001
        print(f'Failed esearch: {e}')
        uids = []
    print(f'Found {len(uids)} nuccore UIDs')
    acc_meta = esummary_accessions(uids) if uids else []
    print(f'Retrieved metadata for {len(acc_meta)} records')
    # Filter length + title includes complete genome
    filtered = [t for t in acc_meta if 29000 <= t[1] <= 31000 and 'complete genome' in t[2].lower()]
    print(f'{len(filtered)} records pass length & title filter (29k-31k & complete genome)')
    # Deduplicate by accession
    unique: dict[str, tuple[str, int, str]] = {}
    for acc, ln, title in filtered:
        unique.setdefault(acc, (acc, ln, title))
    acc_list = sorted(unique.keys())[:MAX_ACCESSIONS]
    if not acc_list:
        # Static minimal fallback: reference isolate Wuhan-Hu-1 (may update version over time)
        acc_list = ['NC_045512.2']
        print('Using static fallback accession list (reference genome)')
    print(f'Selected {len(acc_list)} accession versions (cap {MAX_ACCESSIONS})')

    # Write selection for reproducibility
    list_path = HERE / 'selected_accessions.txt'
    static_path = HERE / 'static_accessions.txt'
    for path in (list_path, static_path):
        with path.open('w') as outf:
            for acc in acc_list:
                outf.write(acc + '\n')
    print(f'Wrote accession list to {list_path} and {static_path}')

    kept = 0
    for i, acc in enumerate(acc_list, 1):
        try:
            fasta_txt = fetch_fasta(acc)
            header = fasta_txt.splitlines()[0] if fasta_txt else ''
            if 'complete genome' not in header.lower():
                print(f'Skip {acc}: header lacks complete genome')
            else:
                # Evaluate N content (exclude header lines starting with '>')
                seq_lines = [l.strip() for l in fasta_txt.splitlines() if l and not l.startswith('>')]
                seq = ''.join(seq_lines)
                if not seq:
                    print(f'Skip {acc}: empty sequence')
                else:
                    n_count = seq.upper().count('N')
                    frac = n_count / len(seq)
                    if frac >= MAX_N_FRACTION:
                        print(f'Skip {acc}: N fraction {frac:.2%} >= {MAX_N_FRACTION:.0%}')
                    else:
                        out_path = HERE / f'{acc}.fasta.gz'
                        with gzip.open(out_path, 'wt') as gz:
                            gz.write(fasta_txt if fasta_txt.startswith('>') else f'>{acc}\n{fasta_txt}')
                        kept += 1
        except Exception as e:  # noqa: BLE001
            print(f'Failed {acc}: {e}')
        sleep(0.34)
        if i % 5 == 0 or i == len(acc_list):
            print(f'Processed {i}/{len(acc_list)} (kept {kept})')
    print(f'Finished: wrote {kept} compressed FASTA files.')


if __name__ == '__main__':  # pragma: no cover
    main()
