#!/usr/bin/env python2.7
"""
    Extract bit scores and sequence lengths from a file created with blast
    -outfmt "7 std qlen slen". Take only the highest score per pair of
    sequences.

    extract_blast_fmt7_bitscore.py <input file>
"""

import sys

hits = {}
with open(sys.argv[1], 'rb') as ih:
    for line in ih:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        seq0, seq1 = fields[:2]
        bitscore = float(fields[11])
        ql = fields[12]
        sl = fields[13]
        key = (seq0, seq1)
        if key not in hits or hits[key] < bitscore:
            hits[key] = (bitscore, ql, sl)

for (s0, s1), (bt, ql, sl) in sorted(hits.iteritems()):
    sys.stdout.write(s0 + '\t' + s1 + '\t' + str(bt) + '\t' + ql + '\t' + sl + '\n')
