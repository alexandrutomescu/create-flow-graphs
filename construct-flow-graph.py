from math import ceil, floor
import sys
import networkx as nx
from graphviz import Digraph
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

rng = default_rng()
genomeAbundances = [ceil(x*1000) for x in rng.lognormal(mean=-4, sigma=4, size = len(genomeFiles))]

def nextKmer(kmer, newBase):
    return kmer[-(len(kmer) - 1):] + newBase

def prevKmer(kmer, newBase):
    return newBase + kmer[:len(kmer) - 1]

def out_neighbors(graph, kmer):
    n = [] # The list of out-neighbors
    for base in dnaBases:
        nKmer = nextKmer(kmer, base)
        if nKmer in graph:
            n.append(nKmer)
            
    return n

def get_genome(filePath):
    genome = ''
    fastqLines = open(filePath, 'r').readlines()

    for lineIndex, line in enumerate(fastqLines):
        if lineIndex == 0:
            continue
        genome += line.strip().upper()
        
    return genome

def augment_dbGraph_from_string(string, order, abundance):
    global dbGraph

    for i in range(len(string) - (order-1)):
        kmer = string[i:i+order]
        dbGraph[kmer] = dbGraph.get(kmer, 0) + abundance

def convert_to_networkx_graph(graph):
    global nextId, s, t
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

    for v in G.nodes():
        if v == s or v == t:
            continue
        in_flow = 0
        for u in G.predecessors(v):
            in_flow += G[u][v]["weight"]
        out_flow = 0
        for w in G.successors(v):
            out_flow += G[v][w]["weight"]

        if in_flow < out_flow:
            G.add_edge(s, v, weight=out_flow - in_flow)
        if in_flow > out_flow:
            G.add_edge(v, t, weight=in_flow - out_flow)
    
    return G

def compact_unary_nodes(G):
    
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
            G.add_edge(u, w, weight=f)

def check_flow_conservation(G):
    
    global s,t
    for node in G.nodes():
        if node not in [s,t]:
            in_flow = sum([G[u][node]["weight"] for u in G.predecessors(node)])
            out_flow = sum([G[node][w]["weight"] for w in G.successors(node)])
            if in_flow != out_flow:
                return False

    return True

def network2dot(nxDigraph):
    dot = Digraph()
    dot.graph_attr['rankdir'] = 'LR' # Display the graph in landscape mode
    dot.node_attr['shape'] = 'rectangle' # Rectangle nodes
    
    for (u,v) in nxDigraph.edges():
        att = nxDigraph[u][v]
        dot.edge(str(u),str(v),label=f'{att["weight"]}')

    return dot

def count_simple_cycles(G, bound):
    count = 0
    for cycle in nx.simple_cycles(G):
        count += 1
        if count == bound:
            break
    
    return count

def write_to_catfish_format(G, filename):
    global genomeFiles
    global genomeAbundances

    f = open(filename,"w")

    f.write('#genomes: ')
    for genomeFile in genomeFiles:
        f.write(f'{genomeFile} ')
    f.write('\n#abundances: ')
    for genomeAbundance in genomeAbundances:
        f.write(f'{genomeAbundance} ')
    f.write('\n')

    f.write(f'{G.number_of_edges()}\n')
    for u,v,a in G.edges(data=True):
        f.write(f'{u} {v} {a["weight"]}\n')
    f.close()


k = 15
gt = 5

genomes = []
range_increment = 2000
min_length = sys.maxsize

for genomeFile in genomeFiles[:gt]:
    genome = get_genome('ecoli/' + genomeFile)
    min_length = min(min_length, len(genome))
    genomes.append(genome)

for range_start in range(0,min_length,range_increment):
    dbGraph = dict() # edges and their abundances
    kmer2id = dict()
    nextId = 2
    s = 0
    t = 1
    for index, genome in enumerate(genomes[:gt]):
        augment_dbGraph_from_string(genome[range_start:range_start+range_increment], k, genomeAbundances[index])

    dbGraph_nx = convert_to_networkx_graph(dbGraph)
    print("Created the graph")
    compact_unary_nodes(dbGraph_nx)
    assert(check_flow_conservation(dbGraph_nx))

    print("Compacted the graph")
    print(f'Number of nodes: {dbGraph_nx.number_of_nodes()}')
    print(f'Number of edges: {dbGraph_nx.number_of_edges()}')
    print(f'Is acyclic: {nx.is_directed_acyclic_graph(dbGraph_nx)}')
    n_cycles = count_simple_cycles(dbGraph_nx, 1000)
    print(f'Truncated number of cycles: {n_cycles}')

    if n_cycles < 50:
        continue
    filename = f'graphs/gt{gt}.kmer{k}.({range_start}.{range_start+range_increment}).V{dbGraph_nx.number_of_nodes()}.E{dbGraph_nx.number_of_edges()}.cyc{n_cycles}.graph'
    write_to_catfish_format(dbGraph_nx, filename)
    # dbGraph_dot = network2dot(dbGraph_nx)
    # dbGraph_dot.view()