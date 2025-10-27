from math import ceil, floor, sqrt, exp
from re import L
import sys
import os
import shlex
import networkx as nx
from graphviz import Digraph
import argparse
from numpy.random import default_rng
try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None
from datasets import load_dataset
import gzip

# Will be assigned based on selected dataset (ecoli or labmix) by listing directory contents.
genomeFiles = []


def get_genome(filePath):
    """Return concatenated uppercase sequence(s) from (possibly gzipped) FASTA file using Biopython if available."""
    is_gz = filePath.endswith('.gz')
    open_func = gzip.open if is_gz else open
    mode = 'rt'  # text mode for gzip as well
    if SeqIO is not None:
        seqs = []
        with open_func(filePath, mode) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seqs.append(str(record.seq).upper())
        return ''.join(seqs)
    # Fallback manual parser (simple, ignores description lines beginning with '>').
    genome = []
    with open_func(filePath, mode) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                continue
            genome.append(line.strip().upper())
    return ''.join(genome)

def augment_dbGraph_from_string(string, index, order, abundance):
    global dbGraph, kmer_paths

    for i in range(len(string) - (order-1)):
        kmer = string[i:i+order]
        dbGraph[kmer] = dbGraph.get(kmer, 0) + abundance
        kmer1 = kmer[:-1]
        kmer2 = kmer[1:]
        if len(kmer_paths[index]) == 0:
            kmer_paths[index].append(kmer1)
        kmer_paths[index].append(kmer2)


def construct_subpaths(genome, index, order):
    global read_starts, subpaths, kmer2id

    subpaths[index] = []
    for start in read_starts[index]:
        new_path = []
        for i in range(min(len(genome)-start-1-order, args.readlength)):
            # nodes are k-1 mers
            kmer = genome[start+i:start+i+order-1]
            nodeid = kmer2id[kmer]
            new_path.append(nodeid)
        subpaths[index].append(new_path)

def convert_to_networkx_graph(graph):
    global nextId, s, t, kmer2id, paths, kmer_paths, genomeAbundances
    G = nx.DiGraph()
    
    for key in graph.keys():
        kmer1 = key[:-1]
        kmer2 = key[1:]
        for kmer in [kmer1, kmer2]:
            if kmer not in kmer2id:
                kmer2id[kmer] = nextId
                nextId += 1
        G.add_edge(kmer2id[kmer1], kmer2id[kmer2], weight=graph[key])
    G.add_nodes_from([s, t])

    # construct paths
    for path_index, kmer_path in enumerate(kmer_paths):
        paths[path_index] = [s]
        for kmer in kmer_path:
            kmerid = kmer2id[kmer]
            paths[path_index].append(kmerid)
        paths[path_index].append(t)

    for index, path in enumerate(paths):
        if G.has_edge(s,path[1]):
            G[s][path[1]]["weight"] += genomeAbundances[index]
        else:
            G.add_edge(s, path[1], weight=genomeAbundances[index])
        
        if G.has_edge(path[-2],t):
            G[path[-2]][t]["weight"] += genomeAbundances[index]
        else:
            G.add_edge(path[-2], t, weight=genomeAbundances[index])

    return G

def compact_unary_nodes(G):
    global removedNodes, subpaths

    unaryNodes = []
    for v in G.nodes():
        if G.in_degree(v) == 1 and G.out_degree(v) == 1:
            unaryNodes.append(v)
    for node in unaryNodes:
        u = list(G.predecessors(node))[0]
        w = list(G.successors(node))[0]
        if not G.has_edge(u,w):
            f = G[u][node]["weight"]
            G.remove_node(node)
            removedNodes.add(node)
            G.add_edge(u, w, weight=f)
            for gen_index in range(len(subpaths)):
                for subpath_index in range(len(subpaths[gen_index])):
                    if len(subpaths[gen_index][subpath_index]) > 0:
                        if node == subpaths[gen_index][subpath_index][0]:
                            subpaths[gen_index][subpath_index][0] = u
                        if node == subpaths[gen_index][subpath_index][-1]:
                           subpaths[gen_index][subpath_index][-1] = w

def satisfies_flow_conservation(G, tol=1e-6, verbose=False):
    """Return True if every internal node satisfies flow conservation within tolerance.

    tol: absolute tolerance for comparing floating sums (needed after introducing float abundances).
    verbose: if True, print first few violations for debugging.
    """
    global s, t
    violations = 0
    for node in G.nodes():
        if node in (s, t):
            continue
        in_flow = sum(G[u][node]["weight"] for u in G.predecessors(node))
        out_flow = sum(G[node][w]["weight"] for w in G.successors(node))
        if abs(in_flow - out_flow) > tol:
            violations += 1
            if verbose and violations <= 5:
                print(f"Flow violation at node {node}: in={in_flow} out={out_flow} diff={in_flow-out_flow}")
            if violations > 5 and verbose:
                print("... further violations suppressed ...")
                break
    return violations == 0

def check_subpaths(G):
    global s, t, subpaths
    for gen_index in range(len(subpaths)):
        for subpath_index in range(len(subpaths[gen_index])):
            subpath = subpaths[gen_index][subpath_index]
            # check that every edge in subpath (i.e. pairs of adjacent list entries, i.e. nodes) is an edge of G
            for i in range(len(subpath)-1):
                if not G.has_edge(subpath[i], subpath[i+1]):
                    print(f"Invalid subpath edge {(subpath[i], subpath[i+1])} for genome {gen_index}: {subpath}")
                    return False
    return True

def network2dot(nxDigraph, style="default", left_index=None, right_index=None):
    # Compute nodes that are reachable from left_index and can reach right_index
    nodes_to_include = set(nxDigraph.nodes())
    
    if left_index is not None:
        # Get nodes reachable from left_index (forward reachability)
        if left_index in nxDigraph:
            reachable_from_left = nx.descendants(nxDigraph, left_index)
            reachable_from_left.add(left_index)  # include the source node itself
            nodes_to_include &= reachable_from_left
        else:
            # If left_index doesn't exist in graph, no nodes pass the filter
            nodes_to_include = set()
    
    if right_index is not None:
        # Get nodes that can reach right_index (backward reachability)
        if right_index in nxDigraph:
            can_reach_right = nx.ancestors(nxDigraph, right_index)
            can_reach_right.add(right_index)  # include the target node itself
            nodes_to_include &= can_reach_right
        else:
            # If right_index doesn't exist in graph, no nodes pass the filter
            nodes_to_include = set()
    
    if style == "default":
        dot = Digraph(format='pdf')
        dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
        dot.node_attr['shape'] = 'rectangle' # Rectangle nodes
        
        for (u,v) in nxDigraph.edges():
            # Only include edges where both nodes are in the filtered set
            if u in nodes_to_include and v in nodes_to_include:
                att = nxDigraph[u][v]
                dot.edge(str(u),str(v),label=f'{att["weight"]}')
    elif style == "RCF2026":
        dot = Digraph(format='pdf')
        dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
        dot.node_attr['shape'] = 'point'
        dot.node_attr['label'] = '' 
        dot.node_attr['width'] = '0.1' 
        # dot.edge_attr['fontsize'] = '20pt'
        
        for (u,v) in nxDigraph.edges():
            # Only include edges where both nodes are in the filtered set
            if u in nodes_to_include and v in nodes_to_include:
                att = nxDigraph[u][v]
                dot.edge(str(u),str(v))
    else:
        raise ValueError(f"Unknown style {style}")

    return dot

def render_graph_outputs(dot_graph, base_filename, want_pdf=False, want_png=False):
    """Render selected formats and remove intermediate .dot source file.

    Graphviz's render() writes a source file (we currently name with .dot suffix)
    plus the requested output (e.g. .dot.pdf). We don't want to keep the .dot
    source artifacts, so delete them after rendering. If both PDF & PNG are
    requested we reuse the same .dot path for each render, then delete once.
    """
    source_path = base_filename + ".dot"  # we retain existing naming scheme
    generated_any = False
    if want_pdf:
        dot_graph.format = 'pdf'
        dot_graph.render(source_path, cleanup=False)
        generated_any = True
    if want_png:
        dot_graph.format = 'png'
        dot_graph.render(source_path, cleanup=False)
        generated_any = True
    if generated_any:
        # Remove only the .dot source file; keep produced outputs (.dot.pdf/.dot.png)
        try:
            if os.path.exists(source_path):
                os.remove(source_path)
        except OSError:
            pass
    return generated_any

# counts simple cycles up to bound
def count_simple_cycles(G, bound):
    count = 0
    for cycle in nx.simple_cycles(G):
        count += 1
        if count == bound:
            break
    
    return count

def write_to_catfish_format(G, filename):
    global s, t, genomeFiles, genomeAbundances, paths, subpaths, removedNodes, args_str

    f = open(filename + ".graph","w")
    f.write(f'#{os.path.basename(filename)}\n')
    f.write(args_str + '\n')
    f.write(f'#unique source is {s}, unique sink is {t}\n')
    f.write('#genomes: ')
    for genomeFile in genomeFiles[:len(paths)]:
        f.write(f'{genomeFile} ')
    f.write('\n# ground truth paths, in the format \'weight node1 node2 node3 ... \'')
    for index, path in enumerate(paths):
        f.write(f'\n#T {genomeAbundances[index]} ')
        for node in path:
            if node not in removedNodes:
                f.write(f'{node} ')
    if len(subpaths) > 0:
        f.write('\n# subpath constraints, in the format \'node1 node2 node3 ... \'')
        for gen_index in range(len(paths)):
            f.write(f'\n# {gen_index} subpath constraints for genome')
            for index, subpath in enumerate(subpaths[gen_index]):
                f.write(f'\n#S ')
                nodes_to_print = [node for node in subpath if node not in removedNodes]
                # check that there is a graph edge for each pair of adjacent nodes
                if len(nodes_to_print) > 0:
                    for i in range(len(nodes_to_print)-1):
                        if not G.has_edge(nodes_to_print[i], nodes_to_print[i+1]):
                            print(f"Invalid subpath edge {(nodes_to_print[i], nodes_to_print[i+1])} for genome {gen_index}: {subpath}")
                            break
                    for node in nodes_to_print:
                        f.write(f'{node} ')

    f.write(f'\n{G.number_of_edges()}')
    for u,v,a in G.edges(data=True):
        f.write(f'\n{u} {v} {a["weight"]}')
    f.close()

def simulate_read_starts(genome, index):
    global read_starts, args

    nreads = args.nreads

    for _ in range(nreads):
        # choose a random number between 0 and len(genome)-1
        start = rng.integers(0, len(genome))
        read_starts[index].append(start)

parser = argparse.ArgumentParser(
    description="""
    Creates flow de Bruijn graphs from a collection of ecoli genomes
    """,
    formatter_class=argparse.RawTextHelpFormatter
    )
parser.add_argument('-k', '--kmersize', type=int, default=15, help='The kmer size')
parser.add_argument('-w', '--windowsize', type=int, default=2000, help='The length of the genome windows from which to build the graphs. Use 0 for whole genomes.')
parser.add_argument('-a', '--acyclic', action='store_true', help='Keep only acyclic graphs')
parser.add_argument('-c', '--mincycles', type=int, default=0, help='Keep only graphs with at least this many cycles')
parser.add_argument('-d', '--distribution', type=str, default='lognormal-44', help='lognormal-44 or lognormal11')
parser.add_argument('-g', '--ngenomes', type=int, help='The number of ecoli genomes from which to construct the graph', required=True)
parser.add_argument('-r', '--nreads', type=int, default=0, help='The number of reads to simulate per genome', required=False)
parser.add_argument('-l', '--readlength', type=int, help='The length of the simulated reads', required=False)
parser.add_argument('-o', '--outdir', type=str, default='', help='outputdir', required=True)
parser.add_argument('-p', '--pdf', action='store_true', help='Render PDF')
parser.add_argument('--png', action='store_true', help='Render PNG image')
parser.add_argument('-s', '--style', type=str, default='default', help='Graph rendering style passed to network2dot (default: "default")')
parser.add_argument('--poisson-edge-errors', dest='poisson_edge_errors', action='store_true', help='When set, perturb perfect integer edge flows by Poisson sampling (imperfect edge weights).')
parser.add_argument('-D', '--dataset', type=str, default='ecoli', choices=['ecoli','labmix','complex32','medium20','helicobacter-hepaticus','JGI','ebola','measles','sars_cov_2'], help='Dataset: ecoli (default), labmix (fixed), complex32 (abundances from complex32/nanosim.abundances.tsv), medium20 (parsed like complex32), helicobacter-hepaticus (compressed FASTA .gz), JGI (abundances from JGI/nanosim.abundances.tsv), ebola (abundances from ebola/abundances.tsv), measles (abundances from measles/abundances.tsv), sars_cov_2 (no fixed abundances; uses --distribution unless --abundances supplied)')
parser.add_argument('-A', '--abundances', type=str, default=None, help='Comma-separated list of abundances (floats) for the g genomes; mutually exclusive with --distribution and keeps exact float weights (no Poisson perturbation). Length must equal --ngenomes.')
parser.add_argument('--left_index', type=int, default=None, help='Filter nodes: only render nodes with ID >= left_index (inclusive). If not specified, no lower bound filtering.')
parser.add_argument('--right_index', type=int, default=None, help='Filter nodes: only render nodes with ID <= right_index (inclusive). If not specified, no upper bound filtering.')

args = parser.parse_args()

print('--------------------------------------')
print('Running with the following parameters:')
command = (sys.argv[0] + ' ') + ' '.join(shlex.quote(s) for s in sys.argv[1:])
args_str = ('#' + command + '\n#') + '\n#'.join(f'{k}={v}' for k, v in vars(args).items())
print(args_str)
print('--------------------------------------')

if args.acyclic and args.mincycles > 0:
    print("ERROR: you cannot set both --acyclic and --mincycles")
    exit()

# no range validation needed; flag is boolean

# Parse explicit abundances if provided
provided_abundances = None
if args.abundances is not None:
    try:
        provided_abundances = [float(x) for x in args.abundances.split(',') if x.strip() != '']
    except ValueError:
        print('ERROR: Could not parse --abundances (expect comma-separated floats)')
        sys.exit(1)
    if len(provided_abundances) != args.ngenomes:
        print(f'ERROR: --abundances length ({len(provided_abundances)}) does not match --ngenomes ({args.ngenomes})')
        sys.exit(1)
    # Enforce mutual exclusivity: detect if user explicitly passed distribution flag
    argv_flags = set(sys.argv[1:])
    if '-d' in argv_flags or '--distribution' in argv_flags:
        print('ERROR: --distribution cannot be used together with --abundances. Remove -d/--distribution when specifying explicit abundances.')
        sys.exit(1)

k = args.kmersize
outdir = args.outdir.strip().strip("/")

genome_dir, genomeFiles, provided_abundances_loaded, dataset_msgs = load_dataset(args.dataset, args.ngenomes, provided_abundances)
if provided_abundances_loaded is not None:
    provided_abundances = provided_abundances_loaded
for m in dataset_msgs:
    print(m)

if os.path.isdir(outdir):
    print(f"ERROR: {outdir} already exists")
    exit()

# Create nested output directory path (parents) in one call
os.makedirs(outdir, exist_ok=False)

genomes = []
range_increment = args.windowsize
min_length = sys.maxsize
max_length = 0

for gt in range(args.ngenomes,args.ngenomes+1):
    for genomeFile in genomeFiles[:gt]:
        genome = get_genome(f'{genome_dir}/{genomeFile}')
        min_length = min(min_length, len(genome))
        max_length = max(max_length, len(genome))
        genomes.append(genome)

    progress = -1
    if range_increment == 0:
        range_increment = max_length

    for range_start in range(0,min_length,range_increment):
        rng = default_rng()
        if provided_abundances is not None:
            genomeAbundances = provided_abundances  # use floats exactly as provided
        else:
            if args.distribution == 'lognormal-44':
                genomeAbundances = [ceil(x*100) for x in rng.lognormal(mean=-4, sigma=4, size = args.ngenomes)]
            elif args.distribution == 'lognormal11':
                genomeAbundances = [ceil(x*10) for x in rng.lognormal(mean=1, sigma=1, size = args.ngenomes)]
            else:
                print("ERROR: unknown distribution, set either lognormal-44 (default) or lognormal11")

        # progress printing
        old_progress = progress
        progress = int(range_start / min_length * 100)
        if progress > old_progress and progress % 5 == 0:
            print(f'Progress: {progress}%')
        # progress printing

        dbGraph = dict() # edges and their abundances
        kmer2id = dict()
        removedNodes = set() # stores the ids of the removed nodes
        paths = [[] for x in range(gt)] # an element for each path, stores the list of each node id on the path
        kmer_paths = [[] for x in range(gt)] # an element for each path, stores the list of kmer nodes on the path
        read_starts = [[] for x in range(gt)] # a list for each genome, storing the starting positions of the reads
        subpaths = [[] for x in range(gt)] # a list for each genome, storing the nodes of each read, with starting position in read_starts

        s, t, nextId = 0, 1, 2
        for index, genome in enumerate(genomes[:gt]):
            genome_fragment = genome[range_start:range_start+range_increment]
            augment_dbGraph_from_string(genome_fragment, index, k, genomeAbundances[index])

        dbGraph_nx = convert_to_networkx_graph(dbGraph)

        if args.nreads > 0:
            for index, genome in enumerate(genomes[:gt]):
                genome_fragment = genome[range_start:range_start+range_increment]
                simulate_read_starts(genome_fragment, index)
                construct_subpaths(genome_fragment, index, k)

        compact_unary_nodes(dbGraph_nx)
        if not satisfies_flow_conservation(dbGraph_nx, tol=1e-6, verbose=True):
            raise AssertionError('Flow conservation violated (see messages above). Consider investigating abundance inputs or unary compaction.')

        if args.poisson_edge_errors:
            for u, v, data in list(dbGraph_nx.edges(data=True)):
                perfect = data['weight']
                imperfect = int(rng.poisson(perfect))
                if imperfect == 0 and perfect > 0:
                    # Avoid dropping edges with positive perfect flow due to Poisson(λ) = 0 draw
                    imperfect = 1
                data['perfect_weight'] = perfect
                data['weight'] = imperfect

        if args.acyclic:
            if nx.is_directed_acyclic_graph(dbGraph_nx):
                e_tag = 'imp' if args.poisson_edge_errors else 'perf'
                filename = f'{outdir}/gt{gt}.kmer{k}.({range_start}.{range_start+range_increment}).V{dbGraph_nx.number_of_nodes()}.E{dbGraph_nx.number_of_edges()}.{e_tag}.acyc'
                write_to_catfish_format(dbGraph_nx, filename)
                dbGraph_dot = network2dot(dbGraph_nx, style=args.style, left_index=args.left_index, right_index=args.right_index)
                render_graph_outputs(dbGraph_dot, filename, want_pdf=args.pdf, want_png=args.png)
        else:
            n_cycles = 0
            if args.mincycles > 0:
                n_cycles = count_simple_cycles(dbGraph_nx, 100)
            if n_cycles >= args.mincycles:
                e_tag = 'imp' if args.poisson_edge_errors else 'perf'
                filename = f'{outdir}/gt{gt}.kmer{k}.({range_start}.{range_start+range_increment}).V{dbGraph_nx.number_of_nodes()}.E{dbGraph_nx.number_of_edges()}.mincyc{n_cycles}.{e_tag}'
                write_to_catfish_format(dbGraph_nx, filename)
                dbGraph_dot = network2dot(dbGraph_nx, style=args.style, left_index=args.left_index, right_index=args.right_index)
                render_graph_outputs(dbGraph_dot, filename, want_pdf=args.pdf, want_png=args.png)

        # activate these for debugging
    # dbGraph_dot = network2dot(dbGraph_nx, style=args.style, left_index=args.left_index, right_index=args.right_index)
    # dbGraph_dot.view()
        # quit()