# E. coli Dataset (50-assembly subset)

This directory contains 50 Escherichia coli genome assemblies (FASTA) used by the workflow when selecting the dataset key `ecoli`.

Source / Citation:
These assemblies are a subset of the Zenodo dataset:
Jarno N. Alanko. (2022). 3682 E. coli assemblies from NCBI [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6577996 (latest DOI resolving to version) / Version DOI: https://doi.org/10.5281/zenodo.6577997

Description (from source, abridged):
A collection of 3682 E. coli assemblies downloaded from NCBI circa 2020, prepared to replicate the E. coli dataset referenced in the paper "Succinct colored de Bruijn graphs" (Muggli et al., Bioinformatics, 2017. https://doi.org/10.1093/bioinformatics/btx067).

Subset rationale:
This repository includes only 50 assemblies (the first 50 by sorted filename from the full set) to provide a lighter-weight test dataset. No additional modification to individual sequences was performed.

License:
The original Zenodo dataset is released under Creative Commons Attribution 4.0 International (CC BY 4.0). When using this subset, you must give appropriate credit to the original author (Jarno N. Alanko) and provide a link to the license/DOI.

If you use this dataset:
- Cite the Zenodo record (Alanko 2022) above.
- Optionally cite the Muggli et al. paper if using it for methodological context.

Contents:
- `*.fna.gz`: Compressed FASTA genome files (one assembly per file).
- (No abundance file is provided for this dataset.)

Abundances / usage notes:
- The code (`datasets.py`, dataset key `ecoli`) will simply list these genome files.
- You may supply fixed abundances via the `--abundances` CLI option (length must match the number of genomes requested) or allow downstream steps to generate a distribution if supported.
- Filenames are preserved from the original dataset (with `.fna.gz` extension).

Provenance / integrity:
- Files were taken as-is from the original tarball (`coli3682_dataset.tar.gz`) and then a 50-file subset was retained.
- No re-formatting of FASTA headers was performed.

Notes:
- To expand beyond 50 assemblies, download the full dataset from Zenodo using the DOI above and place additional FASTA files into this directory (ensure they retain their original filenames).
- Ensure any provided abundance list corresponds to the (possibly truncated) set of genome files in lexical order if you rely on that ordering.

