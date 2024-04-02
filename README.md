#    Create flow de Bruijn graphs from a collection of ecoli genomes

This script creates de Bruijn graphs from the first `g` genomes in the directory `ecoli`. This directory contains 50 E. coli reference genomes from the dataset [3682 E. coli assemblies from NCBI](https://doi.org/10.5281/zenodo.6577996). 

A de Bruijn graph of order `k` built on a set of strings has a node for each k-mer, i.e. a string of length `k` in one of the strings (duplicated k-mers induce just a single node of the graph). We add an edge between two k-mers having a suffix-prefix overlap of (k-1) chracters if their merged (k+1)-mer appears in one of the strings. Since the genomes have repeated k-mers, this graph will contain cycles.  Moreover, by construction, each of the `g` genomes can be "traced" as a walk from its first k-mer to its last k-mer. We also add a global source *s* connected to all such start k-mers, and a global sink *t* connected from all such ending k-mers. 

To obtain flow values for the edges, we proceed as follows. First, for each of the `g` genomes, we assign a weigth from log-normal distribution of mean -4 and variance 4 (or mean 1 and variance 1, see parameter `--distribution`), as in this [Bioinformatics package](http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html). Then, we initialize the flow values of all edge as *0*, and for each (k+1)-mer in a genome of weight *w*, we add *w* to the flow vlue of the edge corresponding to that (k+1)-mer. Thus, each genome of weight *w* adds flow *w* to the *s*-*t* walk to which it corresponds in the graph. The resulting edge values are still a flow (statisfy flow conservation), and we can continue with the next genome.

We restrict the above procedure to non-overlapping windows of length `--windowsize` (default 2000) inside each of the `g` genomes. You can use the parameter `--acyclic` to keep only acyclic graphs, or the parameter `--mincycles` to keep only graphs containing at least a minimum number of cycles (these are computed with the NetworkX function [simple_cycles](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.cycles.simple_cycles.html)). 
    
## Options:

```
  -h, --help            show this help message and exit
  -k KMERSIZE, --kmersize KMERSIZE [default 15]
                        The kmer size
  -w WINDOWSIZE, --windowsize WINDOWSIZE 
                        The length of the genome windows from which to build the graphs [default 2000]
  -a, --acyclic
                        Keep only acyclic graphs
  -c MINCYCLES, --mincycles MINCYCLES 
                        Keep only graphs with at least this many cycles [default 0]
  -d DISTRIBUTION, --distribution DISTRIBUTION 
                        lognormal-44 or lognormal11 [default lognormal-44]
  -g NGENOMES, --ngenomes NGENOMES
                        The number of ecoli genomes from which to construct the graph
  -o OUTDIR, --outdir OUTDIR
                        outputdir
  -p, ----pdf
                        Render PDF
```

## Example 

```
    python3 construct-flow-graph.py -o graphs-g15-w4000-acyc/ -k 15 -g 15 --acyclic --windowsize 4000 --pdf
```

The above command will produce a directory `graphs-g15-w4000-acyc`, which will contain de Bruijn graphs of order `-k 15` from the first `-g 15` E.coli genomes, in windows of length `--windowsize 2000`, only `--acyclic`, and also render PDFs for these `--pdf`. One such graph is `gt15.kmer15.(4384000.4388000).V22.E42.acyc.graph` which looks like:

![Example E.coli Graph](ecoli/gt15.kmer15.(4384000.4388000).V22.E42.acyc.graph.dot.pdf.png)

## Requirements:

```
    numpy
    networkx
    graphviz
```
