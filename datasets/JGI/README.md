# JGI Dataset

This directory contains the JGI dataset genomes and associated abundance file (`nanosim.abundances.tsv`).

Source / Citation:
The data (or a derivative subset thereof) originates from resources described in the following publication:
[https://doi.org/10.1038/sdata.2016.81](https://doi.org/10.1038/sdata.2016.81

If you use this dataset, please cite the above paper.

Contents:
- `*.fasta(.gz)`: Genome FASTA files.
- `nanosim.abundances.tsv`: Tab-separated file with genome filenames in column 1 and raw abundance values in column 2 (before the ×10 scaling performed by the pipeline).

Notes:
- Filenames must match exactly the first column of `nanosim.abundances.tsv` (with optional `.gz` compression) for the pipeline to load them.
- Abundances are scaled (×10, rounded) inside the code when constructing graphs (see `datasets.py`, dataset key `JGI`).
