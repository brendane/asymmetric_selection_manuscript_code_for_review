#!/usr/bin/env python3
"""
    Given a tab-delimited file constructed by bedtools from SibeliaZ
    gff output and annotations, create a table of LCB presence/absence
    and annotation.

    Input file columns:
    strain.contig, _, _, start in contig, end in contig, _, strand, _,
    LCB name, gff entry for gene (9 columns), extent of overlap
"""

import argparse
import collections
import csv
import urllib.parse

from Bio import SeqIO
import numpy

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('gff')
args = parser.parse_args()

## Gather information
lcbs_length = collections.defaultdict(list)
blocks_per_strain = collections.defaultdict(lambda: collections.defaultdict(set))
annotations_per_strain = collections.defaultdict(lambda: collections.defaultdict(list))
strains = set()
annotated = False
with open(args.gff, 'rt') as handle:
    rdr = csv.reader(filter(lambda line: not line.startswith('#'), handle), delimiter='\t')
    for row in rdr:
        strain = row[0].split('.')[0]
        strains.add(strain)
        lcb = int(row[8].replace('id=', ''))
        blocks_per_strain[lcb][strain].add((row[0], row[3], row[4]))
        lcbs_length[lcb].append(int(row[4]) - int(row[3]) + 1)
        if len(row) > 9:
            annotated = True
            overlap = int(row[18])
            if row[11] not in ['region', 'repeat_region', 'pseudogene']:
                annotations_per_strain[lcb][strain].append((overlap >= (int(row[13]) - int(row[12]) + 1), row[17]))

if annotated:
    print('lcb\tn_strains\tmean_length\tmin_length\tmax_length\tn_loci\t' + '\t'.join(sorted(strains)) +
          '\t' + '\t'.join(s + '__annotation' for s in sorted(strains)))
else:
    print('lcb\tn_strains\tmean_length\tmin_length\tmax_length\tn_loci\t' + '\t'.join(sorted(strains)))

for lcb in sorted(lcbs_length):

    ## LCB length
    lengths = lcbs_length[lcb]
    mean_length, median_length = numpy.mean(lengths), numpy.median(lengths)
    min_length, max_length = min(lengths), max(lengths)

    ## for each strain, get the number of sequences and process the annotations 
    n_strains = len(blocks_per_strain[lcb])
    annot = {}
    n_blocks = {}
    for strain in strains:
        annot[strain] = []
        if strain not in blocks_per_strain[lcb]:
            n_blocks[strain] = 0
        else:
            n_blocks[strain] = len(blocks_per_strain[lcb][strain])
            if annotated:
                a = {}
                child_par = {}
                parents = set()
                for complete, gff_entry in annotations_per_strain[lcb][strain]:
                    if gff_entry == '.':
                        continue
                    else:
                        tags = dict(x.split('=', 1) for x in gff_entry.split(';'))
                    ID = tags['ID']
                    lt = ID
                    parent = ''
                    gene = ''
                    product = ''
                    name = ''
                    if 'locus_tag' in tags:
                        lt = tags['locus_tag']
                    if 'Parent' in tags:
                        parent = tags['Parent']
                        child_par[ID] = parent
                        parents.add(parent)
                    else:
                        child_par[ID] = ''
                    if 'gene' in tags:
                        gene = urllib.parse.unquote(tags['gene'])
                    if 'product' in tags:
                        product = urllib.parse.unquote(tags['product'])
                    if 'Name' in tags:
                        name = urllib.parse.unquote(tags['Name'])
                    a[ID] = [lt, gene, product, name, complete]
                for child, parent in child_par.items():
                    if child in parents: continue
                    lt, gene, name, product, complete = a[child]
                    if parent != '':
                        lt = a[parent][0]
                    else:
                        lt = a[child][0]
                    c = 'partial_overlap'
                    if complete:
                        c = 'completely_contained'
                    annot[strain].append(lt + ':' + gene + ':' + name + ':' + product + ':(' +
                                         c + ')')

    ## Print results for LCB
    if annotated:
        print('\t'.join(map(str,
                            [lcb, n_strains, round(mean_length, 1), min_length, max_length, sum(n_blocks.values())] +
                            [n_blocks[s] for s in sorted(strains)] +
                            [';; '.join(annot[s]) for s in sorted(strains)])))
    else:
        print('\t'.join(map(str,
                            [lcb, n_strains, round(mean_length, 1), min_length, max_length, sum(n_blocks.values())] +
                            [n_blocks[s] for s in sorted(strains)])))
