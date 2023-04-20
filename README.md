#    Create flow de Bruijn graphs from a collection of ecoli genomes
    
## Options:

```
  -h, --help            show this help message and exit
    -k KMERSIZE, --kmersize KMERSIZE
                        The kmer size
  -w WINDOWSIZE, --windowsize WINDOWSIZE
                        The length of the genome windows from which to build the graphs
  -a, --acyclic         Keep only acyclic graphs
  -c MINCYCLES, --mincycles MINCYCLES
                        Keep only graphs with at least this many cycles
  -d DISTRIBUTION, --distribution DISTRIBUTION
                        lognormal-44 or lognormal11
  -g NGENOMES, --ngenomes NGENOMES
                        The number of ecoli genomes from which to construct the graph
  -o OUTDIR, --outdir OUTDIR
                        outputdir
```

## Example 

```
    python3 construct-flow-graph.py -o graphs-g15-w4000-acyc/ -g 15 -a -w 4000
```