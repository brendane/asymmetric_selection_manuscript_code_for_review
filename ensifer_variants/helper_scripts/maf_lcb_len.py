#!/usr/bin/env python3
"""
    Make a tab-delimited table of LCB lengths from a MAF file.

    maf_lcb_len.py <input maf>
"""

import argparse
import collections
import gzip

MafSeq = collections.namedtuple('MafSeq', 'seq_name start aligned_bases strand contig_length seq')

def parse_maf(ih, print_line=False, handle=None):
    record = {'a':{}, 's':{}}
    for line in ih:
        if print_line:
            handle.write(line)
        if line.startswith('a'):
            if len(record['s']) > 0:
                yield record
            record = {'a':{}, 's':{}}
            record['a'] = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
            record['s'] = []
        elif line.startswith('s'):
            fields = line.strip().split()[1:]
            record['s'].append(MafSeq(fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]), fields[5]))
        else:
            continue
    yield record


parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('maf')
args = parser.parse_args()

## Read MAF
open_fun = open
if args.maf.endswith('.gz'): open_fun = gzip.open
with open_fun(args.maf, 'rt') as handle:
    for rec in parse_maf(handle):
        length = len(rec['s'][0].seq)
        print(rec['a']['label'] + '\t' + str(length))
