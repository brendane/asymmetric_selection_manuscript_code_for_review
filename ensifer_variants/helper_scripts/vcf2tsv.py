#!/usr/bin/env python3
"""
    Transform a VCF file into a tab-delimited file with 1s and 0s for
    genotype calls.

    vcf2tsv.py --output <output file> <input file>

    OPTIONS

        --subset        A bed file with positions to include, default is
                        to include everything
        --rs            Include SNP ID
        --min-gt        Minimum genotyping rate (0)
        --impute        Replace missing genotypes with the mean across
                        the entire sample
        --het-miss      Make heterozygous sites missing.
        --qtcat         Make QTCAT style output (don't use --rs); only works
                        with biallelic variants and a haploid genome; will
                        automatically set --het-miss
"""

#==============================================================================#

import argparse
import re

import vcf

#==============================================================================#

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--output')
parser.add_argument('--subset')
parser.add_argument('--min-gt', default=0.0, type=float)
parser.add_argument('--impute', default=False, action='store_true')
parser.add_argument('--rs', default=False, action='store_true')
parser.add_argument('--het-miss', default=False, action='store_true')
parser.add_argument('--qtcat', default=False, action='store_true')
parser.add_argument('input')
args = parser.parse_args()

qtcat = args.qtcat
het_miss = args.het_miss
if qtcat:
    het_miss = True

if args.subset is None:
    positions = None
else:
    positions = set()
    with open(args.subset, 'rt') as handle:
        for line in handle:
            contig, _, pos = line.strip().split()[:3]
            positions.add((contig, pos))

het_count = 0

with open(args.output, 'wt') as out:
    rdr = vcf.Reader(filename=args.input)
    if args.rs:
        if qtcat:
            out.write('names\tchr\tpos\t' + '\t'.join(rdr.samples) + '\n')
        else:
            out.write('contig\tpos\trs\t' + '\t'.join(rdr.samples) + '\n')
    else:
        if qtcat:
            out.write('chr\tpos\t' + '\t'.join(rdr.samples) + '\n')
        else:
            out.write('contig\tpos\t' + '\t'.join(rdr.samples) + '\n')
    for rec in rdr:
        if (positions is not None) and ((rec.CHROM, str(rec.POS)) not in positions):
            continue
        rs = rec.ID
        gts = []
        missed = 0
        sum_gt = 0
        if qtcat:
            # Collapse MNPs down to the first variable site so that
            # QTCAT doesn't get confused
            allele0 = str(rec.REF)
            allele1 = str(rec.ALT[0])
            if len(allele0) > 1:
                for a, b in zip(allele0, allele1):
                    if a != b:
                        allele0, allele1 = a, b
        for s in rec.samples:
            gt = re.split('\||/', s['GT'])
            pl = len(gt)
            if pl > 1 and (gt[0] == '.' or ((gt[0] != gt[1]) and het_miss)):
                g = 'NA'
                missed += 1
                if (gt[0] != gt[1]):
                    het_count += 1
            elif pl == 1 and gt[0] == '.':
                g = 'NA'
                missed += 1
            else:
                g = sum(1 for c in gt if c == '0') / float(pl)
                sum_gt += g
            if qtcat:
                if pl > 2:
                    raise Exception('Cannot produce QTCAT output with multi-allelic sites')
                if g == 'NA':
                  _g = 'N'
                else:
                    _g = [allele1, allele0][int(g)]
                gts.append(_g + _g)
            else:
                gts.append(str(g))
        if 1 - (missed / float(len(gts))) < args.min_gt:
            continue
        if missed > 0 and args.impute:
            mean_gt = sum_gt / float(len(gts) - missed)
            gts = [g if g != 'NA' else mean_gt for g in gts]
        if args.qtcat:
            if args.rs:
                out.write(rs + '\t')
            out.write(rec.CHROM + '\t' + str(rec.POS))
        else:
            out.write(rec.CHROM + '\t' + str(rec.POS))
            if args.rs:
                out.write('\t' + rs)
        out.write('\t' + '\t'.join(str(g) for g in gts) + '\n')

if het_miss:
    print('Found %i heterozygous calls' % het_count)
