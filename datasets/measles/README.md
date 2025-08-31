Measles virus (MeV) complete genome sequences associated with:

Penedos et al. 2015. Assessment of the Utility of Whole Genome Sequencing of Measles Virus in the Surveillance of Imported Cases in the UK, 2010-2012. PLoS ONE 10(11): e0143081. DOI: 10.1371/journal.pone.0143081. PMID: 26569100.

Contents:
- Individual compressed FASTA files for each linked nucleotide record whose header contains "complete genome": <accession>.fasta.gz
- Downloader script: download_measles.py (discovers accessions via PubMed -> nuccore elink for PMID 26569100; filters for complete measles genomes; writes per-accession files)
- Abundance generator: generate_abundances.py (lognormal sampling similar to ebola dataset)

Usage:
  python datasets/measles/download_measles.py
  python datasets/measles/generate_abundances.py --seed 42
