#!/usr/bin/env python3
"""
    Concatenate the records in a FASTA file with 100 Ns in
    between.

    concatenate_reference_fasta.py <infile> <contig order>

    UPDATE 26 Jan 2018: Added more flexibility so that it can work
    with genomes that have lots of contigs. Also, converts sequence to
    uppercase.
"""

#==============================================================================#

import argparse
from sys import stdout
import re

from Bio import SeqIO

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--gap', default=100, type=int)
parser.add_argument('--file', default=False, action='store_true')
parser.add_argument('infile')
parser.add_argument('order', nargs='+')
args = parser.parse_args()

gap = args.gap

infile = args.infile
idx = SeqIO.index(infile, 'fasta')

sq_order = []
if args.file:
    with open(args.order[0], 'rt') as ih:
        for line in ih:
            sq_order.append(line.strip())
else:
    sq_order = args.order

sq_lens = {}
padding = {}
for id, rec in idx.items():
    sq_lens[id] = len(rec)
p = 0
for sq in sq_order:
    padding[sq] = p
    p += sq_lens[sq] + gap

stdout.write('>genome\n')
first = True
for id in sq_order:
    rec = idx[id]
    if first:
        first = False
    else:
        stdout.write('N' * gap)
    stdout.write(str(rec.seq).upper())
stdout.write('\n')
