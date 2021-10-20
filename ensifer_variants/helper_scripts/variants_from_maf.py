#!/usr/bin/env python3
"""
    Get table of allele presence/absence from an MAF file.
"""

import argparse
import collections
import gzip
import itertools
import sys

MafSeq = collections.namedtuple('MafSeq', 'seq_name start aligned_bases strand contig_length seq')

def parse_maf_haplotypes(ih):
    record = {'a':{}, 'haplotypes':collections.defaultdict(set),
              'strain_count':collections.defaultdict(int)}
    for line in ih:
        if line.startswith('a'):
            if len(record['haplotypes']) > 0:
                yield record
            record = {'a':{}, 'haplotypes':collections.defaultdict(set),
                      'strain_count':collections.defaultdict(int)}
            record['a'] = {x.split('=')[0]:x.split('=')[1] for x in line.strip().split()[1:]}
        elif line.startswith('s'):
            fields = line.strip().split()[1:]
            strain = fields[0].split('.')[0]
            record['strain_count'][strain] += 1
            record['haplotypes'][fields[5]].add(strain)
        else:
            continue
    yield record

def get_vars_from_haplotypes(haplotypes):
    h = list(haplotypes.items())
    l = len(h[0][0])
    ret = {}

    ## First find SNPs
    snps = []
    for i in range(l):
        alleles = collections.defaultdict(set)
        for j in range(len(h)):
            b = (h[j][0][i]).upper()
            if b not in {'A', 'T', 'C', 'G'}:
                continue
            alleles[b].add(j)
        if len(alleles) > 1:
            for a, idxs in alleles.items():
                snps.append((str(i) + '-' + a, idxs))
                ret[str(i+1) + '-' + a + '-snp'] = set()
                for j in idxs:
                    for s in h[j][1]:
                        ret[str(i+1) + '-' + a + '-snp'].add(s)

    ## Then find indels
    indels = collections.defaultdict(set)
    for j in range(len(h)):
        d = False
        s = None
        e = None
        for i in range(l):
            b = (h[j][0][i])
            if (b != '-' or i == l-1) and d:
                e = i
                if e == l-1: e +=1
                indels[(s+1, e)].add(j)
                d = False
            elif b == '-' and not d:
                s = i
                d = True

    ## Produce output table
    for se, idxs in sorted(indels.items()):
        drs = '_'.join(map(str, se)) + '-del'
        irs = '_'.join(map(str, se)) + '-ins'
        ret[drs] = set()
        ret[irs] = set()
        done = set()
        for j in idxs:
            done.add(j)
            for s in h[j][1]:
                ret[drs].add(s)
        for j in range(len(h)):
            if j not in done:
                for s in h[j][1]:
                    ret[irs].add(s)

    return ret



parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('maf')
parser.add_argument('strains')
args = parser.parse_args()


strains = []
with open(args.strains, 'rt') as h:
    for line in h:
        strains.append(line.strip())

sys.stdout.write('contig\tpos\trs\t' + '\t'.join(strains) + '\n')
open_fun = open
if args.maf.endswith('.gz'):
    open_fun = gzip.open
with open_fun(args.maf, 'rt') as ih:
    for record in parse_maf_haplotypes(ih):
        var_table = get_vars_from_haplotypes(record['haplotypes'])
        for rs, strain_list in var_table.items():
            sys.stdout.write(record['a']['label'] + '\t' +
                             rs.split('-')[0] + '\t' +
                             record['a']['label'] + '-' + rs)
            for strain in strains:
                if strain in strain_list:
                    sys.stdout.write('\t1')
                else:
                    sys.stdout.write('\t0')
            sys.stdout.write('\n')
        sys.stdout.write(record['a']['label'] + '\t' +
                         '0' + '\t' +
                         record['a']['label'] + '-copynumber')
        for strain in strains:
            sys.stdout.write('\t' + str(record['strain_count'][strain]))
        sys.stdout.write('\n')
