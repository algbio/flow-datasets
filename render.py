#!/usr/bin/env python3
"""Batch render missing Graphviz PDFs for .graph files.

Given a directory under 'graphs', recurse through all subdirectories.
For every file ending in '.graph', if the corresponding '<base>.dot.pdf' does
not exist (where <base> is the path without the trailing '.graph'), create it
using the same logic as in construct-flow-graph.py (edges labeled by weight).

Rendering each graph is limited to 30 seconds to avoid getting stuck on large
instances. Timeouts are reported and skipped.

Exit status is 0 unless an unexpected exception occurs (individual file
failures are logged but don't abort the run).

Usage examples (note: graphs_dir is a positional argument, don't write graphs_dir=...):
    # Render PDFs for any .graph file under graphs/labmix (recursively)
    python render_missing_pdfs.py graphs/labmix

    # Show what would be rendered without running dot
    python render_missing_pdfs.py graphs/labmix --dry-run

    # Re-render even existing PDFs with a longer timeout
    python render_missing_pdfs.py graphs/labmix --force --timeout 120

    # Process a single subdirectory deeper inside graphs
    python render_missing_pdfs.py graphs/ecoli/some_window_dir

If you accidentally invoke it like:
    python render_missing_pdfs.py graphs_dir=graphs/labmix
you will get a "Directory not found" error because the literal path 'graphs_dir=graphs/labmix' probably doesn't exist.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from typing import List, Tuple, Dict, Any
import concurrent.futures as _futures

import networkx as nx
from graphviz import Digraph
from subprocess import TimeoutExpired  # we'll raise this ourselves via subprocess.run(timeout=...)
import subprocess
import shutil


def parse_graph_file(path: str) -> nx.DiGraph:
    """Parse a *.graph file produced by write_to_catfish_format into a DiGraph.

    Format (after comments starting with '#'):
      <E>                (number of edges)
      <u> <v> <weight>   (repeated E times)
    """
    G = nx.DiGraph()
    data_lines: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            data_lines.append(line)
    if not data_lines:
        raise ValueError(f"No graph data lines found in {path}")
    try:
        expected_edges = int(data_lines[0])
    except ValueError as e:
        raise ValueError(f"First non-comment line must be edge count in {path}") from e
    edge_lines = data_lines[1:1+expected_edges]
    if len(edge_lines) != expected_edges:
        raise ValueError(
            f"Edge count mismatch in {path}: header says {expected_edges}, found {len(edge_lines)} lines"
        )
    for ln in edge_lines:
        parts = ln.split()
        if len(parts) != 3:
            raise ValueError(f"Malformed edge line in {path}: '{ln}'")
        u, v, w = parts
        try:
            u_i = int(u); v_i = int(v)
            # We preserve integer if possible else float
            if '.' in w or 'e' in w.lower():
                w_val = float(w)
            else:
                w_val = int(w)
        except ValueError as e:
            raise ValueError(f"Non-numeric value in edge line '{ln}' of {path}") from e
        G.add_edge(u_i, v_i, weight=w_val)
    return G


def network2dot(nxDigraph: nx.DiGraph) -> Digraph:
    """Replicate network2dot from construct-flow-graph.py (edges labeled with weight)."""
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR'
    dot.node_attr['shape'] = 'rectangle'
    for (u, v) in nxDigraph.edges():
        att = nxDigraph[u][v]
        dot.edge(str(u), str(v), label=f'{att["weight"]}')
    return dot


def main():
    parser = argparse.ArgumentParser(
        description="Render missing .dot.pdf files for .graph files under a directory (recursively)."
    )
    parser.add_argument('graphs_dir', help="Directory under 'graphs' to process (can also be a deeper subdirectory).")
    parser.add_argument('--timeout', type=float, default=30.0, help='Per-file render timeout in seconds (default: 30)')
    parser.add_argument('--force', action='store_true', help='Re-render even if .dot.pdf exists')
    parser.add_argument('--dry-run', action='store_true', help='List actions without rendering')
    parser.add_argument('-j', '--jobs', type=int, default=1,
                        help='Parallel jobs (default: 1 = sequential). Use a sensible number (e.g., <= CPU cores).')
    args = parser.parse_args()

    root = args.graphs_dir.rstrip('/')
    if not os.path.isdir(root):
        print(f"ERROR: Directory not found: {root}", file=sys.stderr)
        sys.exit(1)

    to_process = []
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if not fname.endswith('.graph'):
                continue
            graph_path = os.path.join(dirpath, fname)
            base = graph_path[:-6]  # strip .graph
            pdf_path = base + '.dot.pdf'
            if args.force or not os.path.exists(pdf_path):
                to_process.append((graph_path, pdf_path))

    # Sort alphabetically by relative path for deterministic processing order
    to_process.sort(key=lambda pair: os.path.relpath(pair[0], root))

    if not to_process:
        print('Nothing to do (no missing PDFs).')
        return

    print(f"Found {len(to_process)} .graph files needing render.")
    successes = 0
    failures = 0
    timed_out = 0

    # Ensure 'dot' executable is available
    dot_exe = shutil.which('dot')
    if dot_exe is None:
        print("ERROR: 'dot' executable not found in PATH. Install Graphviz.")
        sys.exit(1)

    def _render_one(task: Tuple[str, str]) -> Tuple[str, Dict[str, Any]]:
        """Worker function: returns (graph_path, result_info)."""
        graph_path, pdf_path = task
        base = graph_path[:-6]
        dot_file = base + '.dot'
        start = time.time()
        try:
            G = parse_graph_file(graph_path)
            dot = network2dot(G)
            with open(dot_file, 'w') as df:
                df.write(dot.source)
            cmd = [dot_exe, '-Tpdf', dot_file, '-o', pdf_path]
            subprocess.run(cmd, check=True, timeout=args.timeout, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            elapsed = time.time() - start
            if os.path.isfile(pdf_path):
                return graph_path, {"status": "ok", "elapsed": elapsed}
            return graph_path, {"status": "fail", "message": "PDF not created"}
        except TimeoutExpired:
            # Clean up
            try:
                if os.path.isfile(dot_file):
                    os.remove(dot_file)
            except Exception:
                pass
            try:
                if os.path.isfile(pdf_path) and os.path.getsize(pdf_path) == 0:
                    os.remove(pdf_path)
            except Exception:
                pass
            return graph_path, {"status": "timeout"}
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.decode('utf-8', errors='ignore') if e.stderr else ''
            return graph_path, {"status": "fail", "message": f"dot exit {e.returncode}: {stderr.splitlines()[-1] if stderr else ''}"}
        except Exception as e:  # noqa
            return graph_path, {"status": "fail", "message": str(e)}
        finally:
            try:
                if os.path.isfile(dot_file):
                    os.remove(dot_file)
            except Exception:
                pass

    if args.dry_run:
        for idx, (graph_path, pdf_path) in enumerate(to_process, start=1):
            rel = os.path.relpath(graph_path, root)
            print(f"[{idx}/{len(to_process)}] Would render {rel} -> {os.path.basename(pdf_path)}")
        print('Dry run complete.')
        return

    # Sequential fallback preserves original deterministic order
    if args.jobs <= 1:
        for idx, task in enumerate(to_process, start=1):
            graph_path, pdf_path = task
            rel = os.path.relpath(graph_path, root)
            print(f"[{idx}/{len(to_process)}] Rendering {rel} -> {os.path.basename(pdf_path)}")
            _, info = _render_one(task)
            status = info["status"]
            if status == 'ok':
                print(f"  OK ({info['elapsed']:.2f}s)")
                successes += 1
            elif status == 'timeout':
                print(f"  TIMEOUT after {args.timeout}s; skipped")
                timed_out += 1
            else:
                print(f"  FAILED: {info.get('message','unknown error')}")
                failures += 1
    else:
        jobs = max(1, args.jobs)
        print(f"Running with {jobs} parallel jobs...")
        # Submit all tasks first (sorted order) so indices known
        index_map = {gp: i for i, (gp, _) in enumerate(to_process, start=1)}
        rel_cache = {gp: os.path.relpath(gp, root) for gp, _ in to_process}
        with _futures.ThreadPoolExecutor(max_workers=jobs) as ex:
            future_map = {ex.submit(_render_one, task): task for task in to_process}
            completed = 0
            try:
                for fut in _futures.as_completed(future_map):
                    graph_path, result = fut.result()
                    completed += 1
                    idx = index_map[graph_path]
                    pdf_path = future_map[fut][1]
                    print(f"[{idx}/{len(to_process)}] {rel_cache[graph_path]} -> {os.path.basename(pdf_path)}")
                    status = result['status']
                    if status == 'ok':
                        print(f"  OK ({result['elapsed']:.2f}s)")
                        successes += 1
                    elif status == 'timeout':
                        print(f"  TIMEOUT after {args.timeout}s; skipped")
                        timed_out += 1
                    else:
                        print(f"  FAILED: {result.get('message','unknown error')}")
                        failures += 1
            except KeyboardInterrupt:
                print("\nKeyboardInterrupt received; cancelling remaining tasks...")
                for f in future_map:
                    f.cancel()
                # Drain cancelled futures
                _ = [f for f in future_map if f.cancelled()]
                print("Cancelled. Partial summary below.")


    print('\nSummary:')
    print(f"  Successful renders: {successes}")
    print(f"  Timeouts:           {timed_out}")
    print(f"  Failures:           {failures}")
    if failures or timed_out:
        sys.exit(2)


if __name__ == '__main__':
    main()
