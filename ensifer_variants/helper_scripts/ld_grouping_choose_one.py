#!/usr/bin/env python2
"""
    Choose one variant to represent each group from a file with LD
    groupings. Chooses the first one in the file, so the way the file
    is sorted will affect the choice.

    ld_grouping_choose_one.py --output <output file> <input file>

    The output is in bed format, with columns five and six holding the
    group number and all the variants associated with the group,
    respectively.
"""

#==============================================================================#

import argparse
import collections
import csv

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('input')
args = parser.parse_args()

choices = {}
groupings = collections.defaultdict(list)
with open(args.input, 'rb') as inhandle:
    rdr = csv.DictReader(inhandle, delimiter='\t')
    use_rs = ('rs' in rdr.fieldnames)
    for row in rdr:
        if use_rs:
            key = (row['rs'], row['chrom'], int(row['pos']))
        else:
            key = (row['chrom'], int(row['pos']))
        grp = int(row['group'])
        if grp not in choices:
            choices[grp] = key
        groupings[grp].append(row)

with open(args.output, 'wb') as out:
    for g, k in sorted(choices.iteritems()):
        if use_rs:
            ky = str(k[0])
            c = k[1]
            p = k[2]
        else:
            ky = str(k[0]) + '-' + str(k[1])
            c = k[0]
            p = k[1]
        out.write(str(c) + '\t' + str(p-1) + '\t' + str(p) + '\t' +
                  str(ky) + '\t' + str(g) + '\t')
        if use_rs:
            out.write(','.join(x['rs'] for x in groupings[g]) + '\n')
        else:
            out.write(','.join(x['chrom'] + '-' + str(x['pos'])
                               for x in groupings[g]) + '\n')

