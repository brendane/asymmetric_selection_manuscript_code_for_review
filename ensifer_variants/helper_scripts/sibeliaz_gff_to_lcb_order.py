#!/usr/bin/env python3
"""
    Given a GFF output file from SibeliaZ, produce an ordered list of
    LCBs.

    sibeliaz_gff_to_lcb_order.py --output <output directory> 
        [--lcbs <lcb list file>] <input> <strain list file>
"""

import argparse
import csv
import gzip
import sys

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--lcbs')
parser.add_argument('gff')
parser.add_argument('strains')
args = parser.parse_args()

lcbs_to_keep = set()
if args.lcbs is not None:
    open_fun = open
    if args.lcbs.endswith('.gz'): open_fun = gzip.open
    with open_fun(args.lcbs, 'rt') as h:
        for line in h:
            lcbs_to_keep.add(line.strip())

strains = set()
if args.strains is not None:
    with open(args.strains, 'rt') as h:
        for line in h:
            strains.add(line.strip())

data = []
with open(args.gff, 'rt') as h:
    rdr = csv.reader(h, delimiter='\t')
    for row in rdr:
        if row[0].startswith('#'): continue
        data.append(row)

for strain in strains:

    lcb_list = []
    for row in data:
        if row[0].split('.')[0] == strain:
            ctg = row[0].split('.', 1)[1]
            pos = int(row[3])
            strand = row[6]
            lcb = row[8].replace('id=', '')
            if args.lcbs is None or lcb in lcbs_to_keep:
                lcb_list.append((ctg, pos, strand, lcb))
    lcb_list.sort()

    with open(args.output + '/' + strain + '.tsv', 'wt') as out:
        prev_contig = None
        for c, p, s, l in lcb_list:
            if prev_contig is None or prev_contig != c:
                out.write('CONTIGEND\t0\n')
                prev_contig = c
            out.write(l + '\t' + s + '\n')
        out.write('CONTIGEND\t0\n')
