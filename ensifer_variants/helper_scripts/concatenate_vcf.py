#!/usr/bin/env python3
"""
    Concatenate the matching samples from any number of VCF or BCF files.

    Outputs uncompressed vcf.

    concatenate_vcf.py input files...
"""

import argparse
import sys

import pysam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('infiles', nargs='+')
args = parser.parse_args()

vfs = []
samples = set()
for i, fname in enumerate(args.infiles):
    vfs.append(pysam.VariantFile(fname))
    if i == 0:
        samples.update(vfs[0].header.samples)
    else:
        samples.intersection_update(set(vfs[-1].header.samples))
samples = list(samples)
samples.sort()
for vf in vfs:
    vf.subset_samples(samples)

## Make output handle
vf_out = pysam.VariantFile(sys.stdout, 'w', header=vfs[0].header)
vf_out.close()
print('#@CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
      '\t'.join(samples))

for vf in vfs:
    for rec in vf:

        info_string = ''
        if len(rec.info) == 0:
            info_string += '.'
        else:
            for k, v in rec.info.items():
                try:
                    info_string += k + '=' + ','.join(str(x) for x in v) + ';'
                except:
                    info_string += k + '=' + str(v)

        ## Not sure this is correct...
        filter_string = ''
        if len(rec.filter) == 0:
            filter_string += '.'
        else:
            try:
                for k, v in rec.filter.items():
                    filter_string += k + '=' + ','.join(str(x) for x in v) + ';'
            except:
                filter_string = '.'

        format_string = ':'.join(rec.format.keys())

        sys.stdout.write(rec.chrom + '\t' + str(rec.pos) + '\t' +
                         rec.id + '\t' + rec.ref + '\t' +
                         ','.join(rec.alts) + '\t' + str(rec.qual) +
                         '\t' + filter_string + '\t' + info_string[:-1] + '\t' +
                         format_string)

        for sample in samples:
            s_ = []
            for k in rec.format.keys():
                try:
                    x = ','.join(str(x) for x in rec.samples[sample].get(k))
                except:
                    x = str(rec.samples[sample].get(k))
                if x == 'None':
                    continue
                s_.append(x)
            if len(s_) == 0: s_ = ['.']
            sys.stdout.write('\t' + ':'.join(s_))

        sys.stdout.write('\n')
