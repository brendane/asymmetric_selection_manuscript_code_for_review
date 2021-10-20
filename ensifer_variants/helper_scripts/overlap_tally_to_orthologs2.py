#!/usr/bin/env python3
"""
    Take the results of overlap_tally_by_strain.py and create an ncol
    formatted graph file for import into igraph.

    overlap_tally_to_orthologs.py <threshold> <input file(s)>

    Requires overlap of at least <threshold> of the larger sequence.
"""

import argparse
import collections
import csv
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('threshold', type=float)
parser.add_argument('inputs', nargs='+')
args = parser.parse_args()

threshold = args.threshold

## Create a graph with all the genes encountered in the input files
## and an edge between all pairs of genes that have more overlap
## than the threshold. Edges are directed.
#g = networkx.DiGraph()
genes = {}
for infile in args.inputs:
    with open(infile, 'rt') as handle:
        rdr = csv.DictReader(handle, delimiter='\t')
        for row in rdr:
            ## Add all genes encountered in the input files
            g1 = row['target_strain'] + ':' + row['target_locus']
            g2 = row['other_strain'] + ':' + row['other_locus']
            if not g1 in genes: genes[g1] = False
            if not g2 in genes: genes[g2] = False
            ## But only add edges if there is enough overlap
            if float(row['target_covered']) >= threshold and row['target_larger'] == '1':
                sys.stdout.write(g1 + '\t' + g2 + '\n')
                genes[g1] = True
                genes[g2] = True


## Print unlinked genes
for g, status in genes.items():
    if not status:
        sys.stdout.write(g + '\t' + g + '\n')

## Next step:
## import igraph
## g = igraph.load('output_of_this_script', format='ncol', directed=False)
## cc = g.clusters(mode=igraph.WEAK)
