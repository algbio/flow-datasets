#!/usr/bin/env python3
"""Download Ebola Makona 2014 complete genomes associated with Gire et al. 2014.

Automatically discovers all nucleotide records linked to the Gire et al. 2014
paper (PMID 25214632) via PubMed -> nuccore elink, filters for headers
containing 'complete genome', and saves each as its own compressed FASTA file
<accession>.fasta.gz. No aggregate multi-FASTA or accession list file is kept.
"""
from pathlib import Path
from time import sleep
from Bio import Entrez
from xml.etree import ElementTree
import ssl
import certifi
import urllib.request
import urllib.error
import sys
import gzip

HERE = Path(__file__).parent
PMID = '25214632'  # Gire et al. 2014

Entrez.email = 'flow-datasets@example.com'  # default email (can be edited if desired)
Entrez.tool = 'create-flow-graphs-downloader'

# Create an SSL context using certifi to avoid local certificate store issues.
_SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())

def _entrez_urlopen(*args, **kwargs):  # type: ignore[no-untyped-def]
    # Inject our context if not provided
    if 'context' not in kwargs:
        kwargs['context'] = _SSL_CONTEXT
    return urllib.request.urlopen(*args, **kwargs)

# Monkey patch Entrez to use our urlopen with proper CA bundle.
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

def fetch_all_accessions_from_pubmed(pmid: str) -> list[str]:
    """Return list of all nucleotide accessions linked to a PubMed ID.

    Uses elink to nuccore then esummary to get accession.version (Caption + Version).
    """
    # Link PubMed -> nuccore IDs
    with Entrez.elink(dbfrom='pubmed', id=pmid, linkname='pubmed_nuccore') as h:
        link_root = ElementTree.parse(h).getroot()
    nuccore_ids = [e.text for e in link_root.findall('.//LinkSetDb/Link/Id')]
    if not nuccore_ids:
        return []
    # Fetch summaries in batches (max 500 ids per call)
    accs: list[str] = []
    for i in range(0, len(nuccore_ids), 200):
        batch = nuccore_ids[i:i+200]
        with Entrez.esummary(db='nuccore', id=','.join(batch)) as h:
            summ_root = ElementTree.parse(h).getroot()
        for doc in summ_root.findall('.//DocSum'):
            # Caption contains accession root, Version has version number
            caption = None
            version = None
            for item in doc.findall('Item'):
                if item.get('Name') == 'Caption':
                    caption = item.text
                elif item.get('Name') == 'Version':
                    version = item.text
            if caption and version:
                accs.append(f"{caption}.{version}")
    # Stable ordering
    accs = sorted(set(accs))
    return accs

def main():
    print('Retrieving accession list via PubMed -> nuccore elink...')
    accessions: list[str] = []
    try:
        accessions = fetch_all_accessions_from_pubmed(PMID)
    except Exception as e:  # noqa: BLE001
        print(f'Warning: elink retrieval failed ({e}). Will attempt static fallback list.')
    if not accessions:
        # Static fallback: expected Makona full set KM233031-KM233129 (inclusive)
        accessions = [f"KM233{n:03d}.1" for n in range(31, 130)]
        print(f'Using static fallback accession list of {len(accessions)} entries (KM233031-KM233129).')
    else:
        print(f'Found {len(accessions)} linked nucleotide records.')
    # Optional sanity message if fewer than expected
    if len(accessions) < 90:
        print('Note: fewer than expected (~99) accessions returned; proceeding with available set.')

    collected: list[str] = []
    for i, acc in enumerate(accessions, 1):
        try:
            seq_txt = fetch_fasta(acc)
            header = seq_txt.splitlines()[0] if seq_txt else ''
            if 'complete genome' not in header.lower():
                print(f'Skipping {acc}: header lacks "complete genome"')
            else:
                collected.append(seq_txt)
                # Write individual compressed FASTA immediately
                indiv_path = HERE / f"{acc}.fasta.gz"
                with gzip.open(indiv_path, 'wt') as gz:
                    gz.write(seq_txt if seq_txt.startswith('>') else f'>{acc}\n{seq_txt}')
        except Exception as e:  # noqa: BLE001
            print(f'Failed {acc}: {e}')
        sleep(0.34)
        if i % 5 == 0 or i == len(accessions):
            print(f'Processed {i}/{len(accessions)} (successful: {len(collected)})')

    print(f'Finished: kept {len(collected)} complete genome sequences (written as <accession>.fasta.gz).')

if __name__ == '__main__':
    main()
