#!/usr/bin/env python3
"""
    Modify a bam file so that it appears as if it was aligned to a
    genome in which the contigs were concatenated and separated by
    100 Ns. Also changes all the chromosome names to 'genome'.

    concatenate_reference_bam.py [--file] [--gap (100)] <contig order>

    Writes to stdout and strips off all the tags.

    UPDATE 26 Jan 2018: Now can read a file with contig order for
    genomes with lots of contigs.
"""

#==============================================================================#

import argparse
import copy
import re
import sys

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--gap', type=int, default=100)
parser.add_argument('--file', default=False, action='store_true')
parser.add_argument('order', nargs='+')
args = parser.parse_args()

gap = args.gap
sq_order = args.order
sq_lens = {}
padding = {}

if args.file:
    sq_order = []
    with open(args.order[0], 'rt') as ih:
        for line in ih:
            sq_order.append(line.strip())

for line in sys.stdin:
    if line.startswith('@'):
        if line.startswith('@SQ'):
            fields = line.strip().split('\t')
            sn = [re.sub('SN:', '', x) for x in fields if x.startswith('SN')][0]
            ln = [re.sub('LN:', '', x) for x in fields if x.startswith('LN')][0]
            sq_lens[sn] = int(ln)
            continue
        else:
            sys.stdout.write(line)
            continue
    else:
        if len(padding) == 0:
            p = 0
            for sq in sq_order:
                padding[sq] = p
                p += sq_lens[sq] + gap
            total_len = sum(sq_lens.values())
            sys.stdout.write('@SQ\tSN:genome\tLN:' + str(total_len) + '\n')

        fields = line.strip().split('\t')
        flag = int(fields[1])
        replicon = fields[2]
        pos = fields[3]
        mreplicon = fields[6]
        mpos = fields[7]
        if mreplicon == '=':
            mreplicon = replicon

        al = (flag & 0x4) == 0
        mal = (flag & 0x8) == 0

        if al:
            new_pos = str(int(pos) + padding[replicon])
            new_replicon = 'genome'
        else:
            new_pos = pos
            new_replicon = 'genome'
        if mal:
           new_mpos = str(int(mpos) + padding[mreplicon])
           new_mreplicon = '='
        else:
            new_mpos = mpos
            new_mreplicon = '='

        new_fields = fields
        new_fields[2] = new_replicon
        new_fields[3] = new_pos
        new_fields[6] = new_mreplicon
        new_fields[7] = new_mpos
        sys.stdout.write('\t'.join(new_fields[:11]) + '\n')
