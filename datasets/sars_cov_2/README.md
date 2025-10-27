SARS-CoV-2 (Severe acute respiratory syndrome coronavirus 2) complete genome sequences.

Contents:
- Individual compressed FASTA files for each selected nucleotide record whose header contains "complete genome": <accession>.fasta.gz
- Downloader script: download_sars_cov_2.py (queries NCBI nuccore for SARS-CoV-2 complete genomes; retains at most 100 distinct accessions; writes per-accession compressed FASTA files).

Selection strategy:
- Search term: ("Severe acute respiratory syndrome coronavirus 2"[Organism]) AND (complete genome) AND ("2020/01/01"[PDAT] : "3000"[PDAT])
- Fetch summaries, collect accession.version, filter length between 29,000 and 31,000 nt (expected SARS-CoV-2 genome size) and presence of "complete genome" in title/definition.
- Deterministic ordering: sort accessions lexicographically and take first 100.

Provenance note: This is a snapshot script-driven selection for reproducible benchmarking, not an exhaustive nor curated epidemiological set.
