# labmix Dataset (5-virus-mix)

This directory contains the labmix (5-virus-mix) dataset genomes and abundance file (`abundances.tsv`).

Source / Citation:
The dataset originates from the study described in the following paper:
http://nar.oxfordjournals.org/content/early/2014/06/27/nar.gku537

Original data repository:
https://github.com/cbg-ethz/5-virus-mix/

If you use this dataset, please cite the above paper and acknowledge the original repository.

Abundances provenance:
The abundances in `abundances.tsv` are taken from Table 1 (PacBio row) of the cited paper, rounded to the nearest integer.

Contents:
- `*.fasta(.gz)`: Genome FASTA files used by the pipeline.
- `abundances.tsv`: Tab-separated file (filename TAB abundance). These abundances are used directly (no scaling) by the pipeline (see `datasets.py`, dataset key `labmix`).

Notes:
- Filenames in the first column of `abundances.tsv` must match the FASTA filenames (optionally with `.gz`).
