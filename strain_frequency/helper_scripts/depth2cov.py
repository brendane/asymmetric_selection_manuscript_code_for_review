#!/usr/bin/env python3
"""
    Convert the output of samtools depth to a coverage per contig.

    Use -aa option to samtools depth for an accurate count.
"""

import sys

contig = None
reads = 0
bases = 0
for line in sys.stdin:
    c, p, r = line.strip().split()
    p, r = int(p), int(r)

    if contig is not None and c != contig:
        print(contig + '\t' + str(float(reads) / bases))
        reads = 0
        bases = 0

    contig = c
    reads += r
    bases += 1

print(contig + '\t' + str(float(reads) / bases))
