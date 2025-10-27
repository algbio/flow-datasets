# complex32 Dataset

This directory contains genomes and abundance file (`nanosim.abundances.tsv`) for the complex32 dataset.

Source / Citation:
The data (or derived subset) is from / described in:
[https://doi.org/10.3389/fmicb.2022.981458](https://doi.org/10.3389/fmicb.2022.981458)

If you use this dataset, please cite the above publication.

Contents:
- `*.fasta(.gz)`: Genome FASTA files.
- `nanosim.abundances.tsv`: Tab-separated file (filename TAB abundance). Abundances are scaled ×10 and rounded inside the pipeline (see `datasets.py`, dataset key `complex32`).

Notes:
- Filenames in column 1 must match the FASTA files (optionally with `.gz` extension) for successful loading.
