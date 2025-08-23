from math import ceil, floor
from re import L
import sys
import os
import pipes
import networkx as nx
from graphviz import Digraph
import argparse
from numpy.random import default_rng

dnaBases = ['A','C','G','T']
genomeFiles = ['GCA_000005845.2_ASM584v2.fna',
    'GCA_000006665.1_ASM666v1.fna',
    'GCA_000007445.1_ASM744v1.fna',
    'GCA_000008865.2_ASM886v2.fna',
    'GCA_000009565.2_ASM956v1.fna',
    'GCA_000010245.1_ASM1024v1.fna',
    'GCA_000010385.1_ASM1038v1.fna',
    'GCA_000010485.1_ASM1048v1.fna',
    'GCA_000010745.1_ASM1074v1.fna',
    'GCA_000010765.1_ASM1076v1.fna',
    'GCA_000013265.1_ASM1326v1.fna',
    'GCA_000013305.1_ASM1330v1.fna',
    'GCA_000014845.1_ASM1484v1.fna',
    'GCA_000017745.1_ASM1774v1.fna',
    'GCA_000017765.1_ASM1776v1.fna',
    'GCA_000017985.1_ASM1798v1.fna',
    'GCA_000019385.1_ASM1938v1.fna',
    'GCA_000019425.1_ASM1942v1.fna',
    'GCA_000019645.1_ASM1964v1.fna',
    'GCA_000021125.1_ASM2112v1.fna',
    'GCA_000022225.1_ASM2222v1.fna',
    'GCA_000022345.1_ASM2234v1.fna',
    'GCA_000022665.2_ASM2266v1.fna',
    'GCA_000023365.1_ASM2336v1.fna',
    'GCA_000023665.1_ASM2366v1.fna',
    'GCA_000025165.1_ASM2516v1.fna',
    'GCA_000025745.1_ASM2574v1.fna',
    'GCA_000026245.1_ASM2624v1.fna',
    'GCA_000026265.1_ASM2626v1.fna',
    'GCA_000026285.2_ASM2628v2.fna',
    'GCA_000026305.1_ASM2630v1.fna',
    'GCA_000026325.2_ASM2632v2.fna',
    'GCA_000026345.1_ASM2634v1.fna',
    'GCA_000026545.1_ASM2654v1.fna',
    'GCA_000027125.1_ASM2712v1.fna',
    'GCA_000091005.1_ASM9100v1.fna',
    'GCA_000146735.1_ASM14673v1.fna',
    'GCA_000147755.2_ASM14775v1.fna',
    'GCA_000147855.3_ASM14785v3.fna',
    'GCA_000148365.1_ASM14836v1.fna',
    'GCA_000148605.1_ASM14860v1.fna',
    'GCA_000155005.1_ASM15500v1.fna',
    'GCA_000155125.1_ASM15512v1.fna',
    'GCA_000157115.2_Escherichia_sp_3_2_53FAA_V2.fna',
    'GCA_000158395.1_ASM15839v1.fna',
    'GCA_000159295.1_ASM15929v1.fna',
    'GCA_000163155.1_ASM16315v1.fna',
    'GCA_000163175.1_ASM16317v1.fna',
    'GCA_000163195.1_ASM16319v1.fna',
    'GCA_000163215.1_ASM16321v1.fna']


def get_genome(filePath):
    genome = ''
    fastqLines = open(filePath, 'r').readlines()

    for lineIndex, line in enumerate(fastqLines):
        if lineIndex == 0:
            continue
        genome += line.strip().upper()
        
    return genome

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
        next_i = 0
        for i in range(min(len(genome)-start-1-order, args.readlength)):
            # nodes are k-1 mers
            kmer = genome[start+i:start+i+order-1]
            nodeid = kmer2id[kmer]
            new_path.append(nodeid)
            next_i = i + 1
        kmer = genome[start+next_i+1:start+next_i+order]
        if len(kmer) == order - 1:
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
                    if node == subpaths[gen_index][subpath_index][0]:
                        subpaths[gen_index][subpath_index][0] = u
                    if node == subpaths[gen_index][subpath_index][-1]:
                        subpaths[gen_index][subpath_index][-1] = w

def satisfies_flow_conservation(G):
    
    global s,t
    for node in G.nodes():
        if node not in [s,t]:
            in_flow = sum([G[u][node]["weight"] for u in G.predecessors(node)])
            out_flow = sum([G[node][w]["weight"] for w in G.successors(node)])
            if in_flow != out_flow:
                return False

    return True

def network2dot(nxDigraph):
    dot = Digraph(format='pdf')
    dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
    dot.node_attr['shape'] = 'rectangle' # Rectangle nodes
    
    for (u,v) in nxDigraph.edges():
        att = nxDigraph[u][v]
        dot.edge(str(u),str(v),label=f'{att["weight"]}')

    return dot

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

    f = open(filename,"w")
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
                if len(nodes_to_print) > 0:
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

args = parser.parse_args()

print('--------------------------------------')
print('Running with the following parameters:')
command = (sys.argv[0] + ' ') + ' '.join( [pipes.quote(s) for s in sys.argv[1:]] )
args_str = ('#' + command + '\n#') + '\n#'.join(f'{k}={v}' for k, v in vars(args).items())
print(args_str)
print('--------------------------------------')

if args.acyclic and args.mincycles > 0:
    print("ERROR: you cannot set both --acyclic and --mincycles")
    exit()

k = args.kmersize
outdir = args.outdir.strip().strip("/")

if os.path.isdir(outdir):
    print(f"ERROR: {outdir} already exists")
    exit()

os.mkdir(outdir)

genomes = []
range_increment = args.windowsize
min_length = sys.maxsize
max_length = 0

for gt in range(args.ngenomes,args.ngenomes+1):
    for genomeFile in genomeFiles[:gt]:
        genome = get_genome('ecoli/' + genomeFile)
        min_length = min(min_length, len(genome))
        max_length = max(max_length, len(genome))
        genomes.append(genome)

    progress = -1
    if range_increment == 0:
        range_increment = max_length

    for range_start in range(0,min_length,range_increment):
        rng = default_rng()
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
        assert(satisfies_flow_conservation(dbGraph_nx))

        if args.acyclic:
            if nx.is_directed_acyclic_graph(dbGraph_nx):
                filename = f'{outdir}/gt{gt}.kmer{k}.({range_start}.{range_start+range_increment}).V{dbGraph_nx.number_of_nodes()}.E{dbGraph_nx.number_of_edges()}.acyc.graph'
                write_to_catfish_format(dbGraph_nx, filename)
                dbGraph_dot = network2dot(dbGraph_nx)
                dbGraph_dot.render(filename + ".dot")
        else:
            n_cycles = 0
            if args.mincycles > 0:
                n_cycles = count_simple_cycles(dbGraph_nx, 100)
            if n_cycles >= args.mincycles:
                filename = f'{outdir}/gt{gt}.kmer{k}.({range_start}.{range_start+range_increment}).V{dbGraph_nx.number_of_nodes()}.E{dbGraph_nx.number_of_edges()}.mincyc{n_cycles}.graph'
                write_to_catfish_format(dbGraph_nx, filename)
                dbGraph_dot = network2dot(dbGraph_nx)
                if args.pdf:
                    dbGraph_dot.render(filename + ".dot")

        # activate these for debugging
        # dbGraph_dot = network2dot(dbGraph_nx)
        # dbGraph_dot.view()
        # quit()