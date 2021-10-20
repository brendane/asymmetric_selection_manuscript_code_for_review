#!/usr/bin/env python3
"""
    Filter a DGRP file.

    filter_dgrp.py --output <output file> --min-gt <minimum gt rate>
        <input file>

    Update 12 Feb 2018: Can subset by position now and by individual.
"""

import argparse
import csv

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--keep')
parser.add_argument('--min-gt', type=float, default=0.)
parser.add_argument('--positions')
parser.add_argument('--order')
parser.add_argument('infile')
args = parser.parse_args()

filter_strains = False
strains = set()
if args.keep is not None:
    filter_strains = True
    with open(args.keep, 'rt') as ih:
        for line in ih:
            strains.add(line.strip().upper())

filter_positions = False
keep = []
if args.positions is not None:
    filter_positions = True
    order = {}
    with open(args.order, 'rt') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            order[row[0]] = int(row[1])
    with open(args.positions, 'rt') as ih:
        rdr = csv.reader(ih, delimiter='\t')
        for row in rdr:
            if row[0] not in order:
                # If some contigs are not included
                # Might be risky
                print('WARNING: %s not included in order file' % row[0])
                continue
            pad = order[row[0]] - 1
            keep.append((int(row[1])+pad, int(row[2])+pad))


with open(args.output, 'wt') as oh:
    with open(args.infile, 'rt') as ih:
        header = ih.readline().strip().split(',')
        strain_cols = [s.upper() for s in header[2:]]
        if filter_strains:
            new_header = ','.join(header[:2] + [s for s in strain_cols if s in strains])
        else:
            new_header = ','.join(header)
        oh.write(new_header + '\n')
        rdr = csv.reader(ih)
        for row in rdr:
            if filter_positions:
                found = False
                p = int(row[0])
                for s, e in keep:
                    if s < p <= e:
                        found = True
                        break
                if not found:
                    continue
            if filter_strains:
                row = row[:2] + [x for x, s in zip(row[2:], strain_cols) if s.upper() in strains]
            n = len(row) - 2
            ngt = sum(1. for b in row[2:] if b.upper() != 'N')
            if ngt / n >= args.min_gt:
                oh.write(','.join(row) + '\n')
