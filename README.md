# flow-datasets

## Repository structure (start here)

Most users only need the ready-made benchmark graph files under `cyclic-graphs/`.

```
cyclic-graphs/
	<dataset-name>/
		perfect-weights/      # Graphs with "perfect" (unperturbed) edge flows
		imperfect-weights/    # Same graphs with stochastic (Poisson / truncated) edge flow perturbations
datasets/                 # Raw genome FASTA files (+ abundance TSVs for some datasets)
construct.py              # Script to (re)generate graphs from datasets
render.py                 # Batch (re)render missing Graphviz PDFs for *.graph files
scripts/                  # Example shell wrappers showing typical construct.py invocations
requirements.txt          # Python dependencies (plus Graphviz system package required for rendering)
README.md                 # You are here
```

### Quick start (just use existing graphs)

If you merely want benchmark flow de Bruijn graphs for evaluation:

1. Pick a dataset directory under `cyclic-graphs/` (e.g. `ecoli`, `labmix`, `complex32`, `JGI`, etc.).
2. Choose `perfect-weights` (deterministic edge flows) or `imperfect-weights` (noisy flows) depending on your experiment.
3. Each `.graph` file encodes one windowed graph; optional `.graph.dot.pdf` (and maybe `.png`) provides a visualization.
4. Parse the `.graph` file to obtain: edge list with weights, ground-truth genome paths (`#T` lines), and optional subpath constraints (`#S` lines).

You only need to run `construct.py` yourself if you want to:
* Change k-mer size, window length, abundance distribution, or noise parameters.
* Generate new graphs for different numbers of genomes or entirely new datasets you supply in `datasets/`.

The remainder of this document explains generation details and parameters.

---

## Installation & Environment

Python 3.10+ recommended.

1. (Recommended) Create a virtual environment.
2. Install Python dependencies:

```bash
pip install -r requirements.txt
```

3. Install Graphviz (needed for rendering):
	* macOS (Homebrew): `brew install graphviz`
	* Ubuntu/Debian: `sudo apt-get install graphviz`
	* Windows (Chocolatey): `choco install graphviz`

4. (Optional) If Biopython fails to install you can still run; a lightweight FASTA parser fallback is used (but Biopython is faster & stricter).

Verify `dot` is on your `PATH`:

```bash
dot -V
```

---

## Output Overview

`construct.py` writes one `.graph` file per processed genome window. Filenames encode metadata, e.g.:

```
gt5.kmer21.(0.2000).V153.E212.mincyc7.e0.5.graph
```

Meaning:

* `gt5` – 5 genomes (ground‑truth paths)
* `kmer21` – k‑mer size k=21
* `(0.2000)` – genomic window [start,end)
* `V153.E212` – 153 nodes, 212 edges after compaction
* `mincyc7` – graph had ≥7 (actually 7) simple cycles (present only if `--mincycles > 0`)
* `acyc` – present instead of `mincyc*` if `--acyclic`
* `e0.5` – imperfect edge flow parameter `--erroreps` (see below)

If PDF / PNG rendering was requested, companion files `.graph.dot.pdf` and/or `.graph.dot.png` are produced (the intermediate `.dot` source is automatically removed in `construct.py`).

### `.graph` File Format

Plain text with comment lines beginning `#`. Key sections:

* `#T <abundance> <node seq...>` – one per ground‑truth path (source `s` to sink `t` inclusive, with any removed unary nodes skipped)
* `#S <node seq...>` – (optional) subpath constraints produced from simulated reads
* Final data block: first a line with the number of edges `E`, followed by `E` lines of `u v weight`.

Edge weights are integer flows. When you supply explicit abundances (`--abundances`) weights are “perfect”. Otherwise, if `--erroreps < 1`, some edges may have an added *imperfect* (truncated Poisson) sample recorded only implicitly via changed `weight` (original `perfect_weight` is added in‑memory only, not in the file).

---

## `construct.py` – Build Flow de Bruijn Graphs

Creates one graph per window (or whole genome) across a chosen number of genomes. Each genome path is embedded; edges are k‑mers (as usual), nodes are (k−1)-mers (after compaction). Unary nodes are contracted while preserving flow.

### CLI Arguments

| Flag | Required | Default | Description |
|------|----------|---------|-------------|
| `-g`, `--ngenomes` | Yes | – | Number of genomes to include (must not exceed dataset size). |
| `-k`, `--kmersize` | No | 15 | k‑mer length `k` used for edges. |
| `-w`, `--windowsize` | No | 2000 | Length of contiguous window taken from each genome. Use `0` for whole genome. Multiple windows are processed sequentially from position 0 up to shortest genome length. |
| `-D`, `--dataset` | No | `ecoli` | Dataset: one of `ecoli`, `labmix`, `complex32`, `medium20`, `helicobacter-hepaticus`, `JGI`, `ebola`, `measles`, `sars_cov_2`. Each may define fixed abundances (see below). |
| `-d`, `--distribution` | No | `lognormal-44` | Abundance distribution when not using dataset‑fixed or explicit abundances. Choices: `lognormal-44` (heavy tailed), `lognormal11` (moderate). Ignored if dataset provides fixed values or if `--abundances` is set. |
| `-A`, `--abundances` | No | `None` | Comma-separated explicit float abundances for exactly `g` genomes. Disallows `--distribution` and disables Poisson perturbation (edge weights remain those values, rounded internally to integers where applied). |
| `-e`, `--erroreps` | No | 1.0 | Imperfect flow sampling epsilon in `[0,1]`. `1.0` = sample full Poisson around each perfect edge flow; `0` = deterministic median (no variance). Values in between truncate Poisson to central interval of width `epsilon` around median before sampling. Ignored if `--abundances` or dataset fixed abundances? (Still applied unless explicit abundances were given; dataset-provided integers are subject to imperfect sampling if `epsilon < 1`). |
| `-a`, `--acyclic` | No | off | Keep only windows whose resulting graph is a DAG (mutually exclusive with `--mincycles`). |
| `-c`, `--mincycles` | No | 0 | Keep only graphs with at least this many (simple) cycles (enumerated up to 100). Mutually exclusive with `--acyclic`. |
| `-r`, `--nreads` | No | 0 | Number of simulated read starts per genome (for subpath constraints). Requires `--readlength`. |
| `-l`, `--readlength` | No | – | Simulated read length in (k−1)-mer node units (actually iterates base positions). Only used if `--nreads > 0`. |
| `-o`, `--outdir` | Yes | – | Output directory (must not already exist). One `.graph` (plus optional renders) per window is created inside. |
| `-p`, `--pdf` | No | off | Emit Graphviz PDF for each graph. |
| `--png` | No | off | Emit Graphviz PNG for each graph. Can be combined with `--pdf`. |

### Dataset Abundance Behavior

Datasets that include an abundances TSV (`labmix`, `complex32`, `medium20`, `JGI`, `ebola`) provide fixed integer (rounded) abundances. In these cases a user-supplied `--abundances` is ignored with an INFO message. For `ecoli` (and any dataset lacking a TSV such as `measles` or `sars_cov_2`) abundances are sampled via `--distribution` unless you pass `--abundances`.

### Imperfect Edge Flows (`--erroreps`)

For each perfect integer edge weight `f`, an “imperfect” weight is sampled from a (possibly truncated) Poisson distribution with mean `f`:

* `epsilon = 1` – sample from full Poisson(f)
* `0 < epsilon < 1` – restrict to central interval of CDF width `epsilon` about the median, then sample
* `epsilon = 0` – deterministic median (equals `f` for integer Poisson median; fallback keeps at least 1 if original was >0)

Set `--erroreps 1` (default) for realistic noise; reduce toward 0 to suppress variation.

### Simulated Read Subpaths

If `--nreads > 0` and `--readlength` provided, random start positions are chosen per genome fragment; each yields a subpath (sequence of node IDs) recorded as `#S ...` lines. These can act as path constraints for downstream reconstruction algorithms.

---

## `render.py` – Batch Produce Missing PDFs

Traverses a directory tree (typically under `graphs/`) and renders a `.dot.pdf` for every `.graph` file lacking one (or all, with `--force`). Rendering is parallelizable and time‑limited per file.

### Arguments

| Positional / Flag | Default | Description |
|-------------------|---------|-------------|
| `graphs_dir` | – | Root directory to scan recursively for `*.graph` files. |
| `--timeout` | 30.0 | Per‑file render timeout in seconds. |
| `--force` | off | Re-render even if a PDF already exists. |
| `--dry-run` | off | List what would be rendered without invoking Graphviz. |
| `-j`, `--jobs` | 1 | Parallel threads (choose ≤ number of CPU cores). |

---

## Provided Convenience Scripts

Scripts in `scripts/` (e.g. `run_ecoli.sh`, `run_labmix.sh`, etc.) are the exact commands that were used to generate the released graphs under `cyclic-graphs/` (with the corresponding parameters and random seeds implicit in each run). They also serve as editable examples if you wish to regenerate or extend the datasets.

If you rerun them today you should obtain structurally comparable graphs; stochastic differences can arise when abundance distributions or Poisson noise (`--erroreps`) are involved. Feel free to duplicate and modify them for custom experiments.

Run one directly (ensure executable bit):

```bash
./scripts/run_ecoli.sh
```

---

## Tips & Troubleshooting

* ERROR about existing output directory: choose a fresh `-o` path; directories are not overwritten.
* Flow conservation assertion failure: indicates a bug or unexpected abundance transform; try re-running with `-e 0` to simplify or inspect earlier log messages.
* Large graphs: rendering may time out; increase `--timeout` or skip rendering during construction and batch render later with `render.py`.
* Missing `dot`: install Graphviz and ensure your shell sees it (`which dot`).

---

Feel free to open issues / PRs for clarifications or enhancements.
