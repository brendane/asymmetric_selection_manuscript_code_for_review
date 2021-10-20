#!/usr/bin/env python3
"""
    Process the GFF output of the SibeliaZ WGA software.
"""

import argparse
import collections
import csv

from Bio import SeqIO

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('gff')
parser.add_argument('fasta')
args = parser.parse_args()

blocks = collections.defaultdict(list)
blocks_count = collections.defaultdict(lambda: collections.defaultdict(int))
blocks_length = collections.defaultdict(list)
strain_set = set()
with open(args.gff, 'rt') as handle:
    rdr = csv.reader(filter(lambda line: not line.startswith('#'), handle), delimiter='\t')
    for row in rdr:
        strain = row[0].split('.')[0]
        strain_set.add(strain)
        block_id = int(row[8].replace('id=', ''))
        blocks_count[block_id][strain] += 1
        blocks_length[block_id].append(int(row[4]) - int(row[3]) + 1)
        blocks[strain].append((row[0], int(row[3])-1, int(row[4])))

i = max(blocks_length.keys())
with open(args.output, 'wt') as oh:
    for strain in sorted(strain_set):
        unaligned = []
        bl = blocks[strain]
        bl.sort()
        contigs = {}
        for rec in SeqIO.parse(args.fasta + '/' + strain + '.fasta', 'fasta'):
            contigs[rec.id] = len(rec)

        contigs_done = set()
        prev_s = 0
        prev_e = 0
        prev_contig = None
        for contig, start, end in bl:
            contigs_done.add(contig)
            if prev_contig == None:
                prev_contig = contig
            if contig != prev_contig:
                if prev_e < contigs[prev_contig]:
                    unaligned.append((prev_contig, prev_e, contigs[prev_contig]))
                #start = 0
                prev_s = 0
                prev_e = 0
                prev_contig = contig
            if start > prev_e:
                unaligned.append((contig, prev_e, start))
            prev_s = start
            prev_e = end
            prev_contig = contig
        if prev_e < contigs[prev_contig]:
            unaligned.append((prev_contig, prev_e, contigs[prev_contig]))
        for contig in set(contigs) - contigs_done:
            unaligned.append((contig, 0, contigs[contig]))
        unaligned.sort()

        for c, s, e in unaligned:
            i += 1
            oh.write('\t'.join([c, 'SibeliaZ', 'singleton_block', str(s+1), str(e), '.', '+', '.',
                                'id=' + str(i)]) + '\n')
