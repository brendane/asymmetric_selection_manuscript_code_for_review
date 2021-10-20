#!/bin/bash
#
# Put PAV and SNP variants from alignments of de novo assembly orthologs
# into a single vcf file.
#
# Notes
# - Ambiguities and gaps are treated as missing data
# - If there are multiple copies of a gene for a strain, differing
#   alleles are treated as ambiguous
# - NCBI annotations for all strains if there was a choice
# - Does not make a DGRP file for HARP, because that analysis makes more
#   sense with a reference genome.
#
# SETTINGS
#
QUEUE=amdsmall
WALLTIME="01:30:00"

PROJDIR="${HOME}/project/select_reseq"
REFERENCE="${PROJDIR}/data/reference/USDA1106/USDA1106.final.fasta"
PAVS="${PROJDIR}/../jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/overlap_orthologs.0.80.tsv"
VCFS="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-04-04/meliloti/vcf"
STRAINS="${PROJDIR}/data/strain/lists/meliloti_sibeliaz_with_freebayes_2020-04-09_ncbi.txt"

OUTDIR="${PROJDIR}/results/combine_variants/meliloti_2021-04-06_denovo"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRIPTDIR="${OUTDIR}/script_copies"

concat="${PROJDIR}/script/bin/concatenate_vcf.py"
vcf2dgrp="${PROJDIR}/script/bin/vcf2dgrp_whole.py"
vcfmissing="${PROJDIR}/script/bin/vcf_missing_as_extra_allele.py"

combine_snp_vcfs="#!/usr/bin/env python3
import collections
import os
import re
import sys
import tempfile

import pysam

strains = {}
with open(sys.argv[1], 'rt') as handle:
    for line in handle:
        x = line.strip()
        strains[x] = x.split('__')[0]

sstrains = sorted(strains)
chroms = []
tf = tempfile.mkstemp(dir='.')[1]
with open(tf, 'wt') as out:
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
              '\t'.join(strains[s] for s in sstrains) + '\n')                 
    n_dups = 0
    for fname in os.listdir(sys.argv[2]):
        with pysam.VariantFile(sys.argv[2] + '/' + fname) as v:
            smap = {}
            for sample in v.header.samples:
                try:
                    smap[sample] = strains[re.sub('(.+)__(.+?)_.+', '\\\\1__\\\\2', sample)]
                except KeyError:
                    continue
            chr = fname.replace('.vcf', '')
            for rec in v:
                dups = False
                if rec.alts[0] == '<*>': continue
                alleles = [rec.ref] + list(rec.alts)
                out.write(chr + '\t' + str(rec.pos) + '\t' +
                          'snp-' + chr + '-' + str(rec.pos) + '\t' +
                          rec.ref + '\t' + ','.join(rec.alts) + '\t' +
                          '.\t.\t.\tGT')
                gts = collections.defaultdict(set)
                for sample in rec.samples:
                    a = rec.samples[sample].alleles
                    if len(a) == 0 or a[0] is None or a[0] == '-' or a[0].upper() == 'N':
                        continue
                    if sample in smap:
                        gts[smap[sample]].add(a[0].upper())
                for s in sstrains:
                    strain = strains[s]
                    if len(gts[strain]) > 1: dups = True
                    if strain in gts and len(gts[strain]) == 1:
                        out.write('\t' + str(alleles.index(list(gts[strain])[0])))
                    else:
                        out.write('\t.')
                out.write('\n')
                if dups:
                    n_dups += 1
            chroms.append((chr, rec.pos))

sys.stdout.write('##fileformat=VCFv4.2\n')
for chr, length in chroms:
    sys.stdout.write('##contig=<ID=%s,length=%i>\n' % (chr, length))
sys.stdout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n')
with open(tf, 'rt') as handle:
    for line in handle:
        sys.stdout.write(line)

os.unlink(tf)
sys.stderr.write('Number of sites with conflicting alleles: %i\n' % n_dups)
"

pavs2vcf="
import csv
import re
import sys
keep = set()
with open(sys.argv[2], 'rt') as handle:
    keep = set(map(lambda x: x.strip(), handle.readlines()))
with open(sys.argv[1], 'rt') as handle:
    rdr = csv.DictReader(handle, delimiter='\t')
    keep = set.intersection(keep, rdr.fieldnames)
    print('##fileformat=VCFv4.2')
    print('##contig=<ID=PAV,length=1000000>\n')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
          '\t'.join(map(lambda s: re.sub('__.+', '', s).replace('_', '-'), keep)))
    for row in rdr:
        o = row['ortholog']
        p = str(int(o.replace('OL', '')) + 1)
        sys.stdout.write('PAV\t' + p + '\tsibeliaz-' + o + '\tA\tP\t.\t.\t.\tGT')
        for strain in keep:
            if row[strain] == '':
                sys.stdout.write('\t0')
            else:
                sys.stdout.write('\t1')
        sys.stdout.write('\n')
"

if [[ "$PBS_ENVIRONMENT" == PBS_BATCH ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.10.2
    module load htslib/1.9
    module load vcflib/1.0.1

    set -euo pipefail

    cd "$OUTDIR"

    cp "$STRAINS" keep.sibeliaz.txt

#    ## Extract SNPs
#    echo "$combine_snp_vcfs" > combine_snp_vcfs.py
#    chmod u+x combine_snp_vcfs.py
#    ./combine_snp_vcfs.py keep.sibeliaz.txt "$VCFS" \
#        | sed 's/_/-/g' \
#        | bcftools norm -Ov -m "-any" \
#        | bcftools annotate -Ov --set-id 'snp-%CHROM-%POS-%FIRST_ALT' \
#        > snps.vcf \
#        || { echo "making SNP vcf file failed"; exit 1; }
#
#    ## Filter and reformat the PAVs
#    python -c "$pavs2vcf" "$PAVS" keep.sibeliaz.txt \
#        | sed 's/_/-/g' \
#        > pavs.vcf \
#        || { echo "processing PAVs failed"; exit 1; }

    ## Merge
    cp "$concat" .
    "./$(basename "$concat")" snps.vcf pavs.vcf \
        | grep -v "^#CHROM" \
        | sed 's/^#@CHROM/#CHROM/g' \
        | bcftools view -Ov --min-ac 1:nonmajor \
        > variants.vcf \
        || { echo "concatenating files failed"; exit 1; }

    ## Compress, index, and tabix; remove info field from bcf file
    bgzip -f -i variants.vcf && tabix -p vcf variants.vcf.gz \
        || { echo "bgzipping and tabixing failed"; exit 1; }
    zcat variants.vcf.gz \
        | vcfkeepinfo - "" \
        | vcfkeepgeno - GT \
        > tmp.vcf \
        || { echo "trimming fields from vcf failed"; exit 1; }
    bgzip -f -i tmp.vcf && tabix -p vcf tmp.vcf.gz \
        || { echo "bgzipping & tabixing tmp failed"; exit 1; }
    bcftools view -Ob tmp.vcf.gz > variants.bcf \
        || { echo "converting to bcf failed"; exit 1; }

#    rm -f pavs.vcf* snps.vcf tmp* snp.tmp.vcf

    echo "RUN TIME = ${SECONDS} ($(($SECONDS/60)))"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$SCRIPTDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --time=$WALLTIME "$sfile"
    echo "$WORKDIR"

fi
