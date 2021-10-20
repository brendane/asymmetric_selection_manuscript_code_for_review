#!/bin/bash
#
# Compare popgen summary stats for candidate and background genes.
#
# Same as the 2021-03-19_libsequence run, but the input data were updated.
# Specifically, fragment and pseudogene sequences were removed more
# thoroughly.
#
# NOTE: The input data for the alignment step that preceeds this code
# (2021-04-04 Pagan aligner run) included a bunch of singleton genes
# that are fragments, which were filtered out before this step. It also
# is missing a few genes from the NCBI PacBio strains, which were
# incorrectly filtered out. However, that will not have any effect
# on the MAG strains here. It is possible that re-running the alignment
# with the missing sequences included would slightly change the positions
# of indels, but this should be quite rare.
#
# UPDATE 23 July 2021: Missing genes added to BetaScan results 
#
# UPDATE 30 July 2021: gwa_n_strains and gwa_n_genes columns fixed so that
# they actually are gwas strains and not all strains.
#
# UPDATE 02 Aug 2021: Added gestimator code to calculate KaKs on medicae/
# meliloti sequences.
#
# SETTINGS
#
QUEUE=amdsmall

PROJDIR="${HOME}/project/popgen2"
ALIGNMENT_DIRECTORY="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-04-04/meliloti"
REPLICONS="${PROJDIR}/../jgi_sequencing/results/replicon_assignments/2021-04-05/meliloti/replicon_assignments.tsv"
ORTHOLOGS="${PROJDIR}/../jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/overlap_orthologs.0.80.tsv"
VCF="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-04-04/meliloti/spanfran_vcf"
VCF2="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-04-04/meliloti/vcf"
BBH="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-06-07/bbh/bbh.tsv"
BBH_ALIGNS="${PROJDIR}/../jgi_sequencing/results/gene_aligns/from_sibeliaz/pagan2/2021-06-07/bbh/aligned-"
BETA="${PROJDIR}/data/tuomas_selection/2021-07-23/betascan_rhizobia.txt"
RAISD="${PROJDIR}/data/tuomas_selection/2021-07-23/raisd_rhizobia.txt"
GWA_STRAINS="${PROJDIR}/../select_reseq/data/strain/lists/spanfran2_meliloti.txt"

OUTDIR="${PROJDIR}/results/popgen_summaries/2021-04-05_libsequence"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYDIR="${OUTDIR}/arrayjobdata"

maf_table="#!/usr/bin/env python3
import csv
import os
import sys

N_STRAINS=165

print('gene\tchrom\tpos\trs\tpav\tnonsyn\tsyn\tspanfran_alt_af\tgwa_alt_af')

gwa_strains = set()
with open('${GWA_STRAINS}', 'rt') as handle:
    for line in handle: gwa_strains.add(line.strip())

ns = set()
sy = set()
with open('nsyn_sites.txt', 'rt') as handle:
    for line in handle:
        ns.add(tuple(line.strip().split('\t')))
with open('syn_sites.txt', 'rt') as handle:
    for line in handle:
        sy.add(tuple(line.strip().split('\t')))

with open('${ORTHOLOGS}', 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        gene = row['ortholog']
        rs = 'sibeliaz-' + gene
        s = set()
        gwa_s = set()
        for st, gs in row.items():
            strain = st.split('__')[0]
            if strain.startswith('MAG') and gs != '':
                s.add(strain)
                if strain in gwa_strains:
                    gwa_s.add(strain)
        print(gene, 'PAV', int(gene.replace('OL',''))+1, rs,
              1, 0, 0, len(s)/N_STRAINS, len(gwa_s)/len(gwa_strains),
              sep='\t')

def div0(x, y):
    if y == 0:
        return float('nan')
    else:
        return x/y

for fname in os.listdir(sys.argv[1]):
    with open(sys.argv[1] + '/' + fname, 'rt') as handle:
        strain_list = []
        for line in handle:
            if line.startswith('##'): continue
            row = line.strip().split('\t')
            if line.startswith('#'):
                for strain in row[9:]:
                    strain_list.append(strain.split('__')[0])
                continue
            chrom, pos = row[:2]
            if row[4] == '.': continue
            for i, alt in enumerate(row[4].split(',')):
                gene = fname.replace('.vcf', '')
                syn = int((gene,pos) in sy)
                nonsyn = int((gene,pos) in ns)
                rs = 'snp-' + gene + '-' + pos + '-' + alt
                counts = 0
                gwa_counts = 0
                n = 0
                gwa_n = 0
                for j, gt in enumerate(row[9:]):
                    if not strain_list[j].startswith('MAG'):
                        continue
                    n += 1
                    if strain_list[j] in gwa_strains:
                        gwa_n += 1
                    if gt == str(i+1):
                        counts += 1
                        if strain_list[j] in gwa_strains:
                            gwa_counts += 1
                if counts > 0:
                    print(gene, chrom, pos, rs, 0, nonsyn, syn,
                          counts / n,
                          div0(gwa_counts, gwa_n), sep='\t')
"

maf="
import collections
import os
import sys
print('gene\tmaf')
for fname in os.listdir(sys.argv[1]):
    maf = []
    with open(sys.argv[1] + '/' + fname, 'rt') as handle:
        for line in handle:
            if line.startswith('#'): continue
            counts = collections.defaultdict(int)
            n = 0
            for gt in line.strip().split('\t')[9:]:
                if gt != '.':
                    counts[gt] += 1
                    n += 1
            if len(counts) > 1:
                maf.append(min(counts.values()) / n)
        maf.sort()
    print(fname, ','.join(map(str, maf)), sep='\t')
"

trim_and_filter="#!/usr/bin/env python3
import os
import re
import sys
from Bio import AlignIO
with open(sys.argv[2], 'rt') as handle:
    strains = set(map(str.strip, handle.readlines()))
for fname in os.listdir(sys.argv[1]):
    if not fname.endswith('.fas'): continue
    a = AlignIO.read(sys.argv[1] + '/' + fname, 'fasta')
    cutat = None
    l = len(a[0])
    N = set('N')
    Ng = set('N-')
    for i in range(-1, -l, -1):
        x = set(a[:,i].upper())
        if x == Ng or x == N:
            cutat = i
        else:
            break
    if cutat is None: exit
    a = a[:,:cutat]
    output = []
    for rec in a:
        if rec.id.split('__', 1)[0] in strains:
            s = str(rec.seq)
            if re.search('N-+$', s) is not None:
                new_s = list(s)
                for i in range(-1, -l, -1):
                    if s[i] == '-': continue
                    elif s[i].upper() != 'N': break
                    else: new_s[i] = '-'
                s = ''.join(new_s)
            output.append((rec.id, s))
    if len(output) == 0: continue
    with open(sys.argv[3] + '/' + fname.replace('.fas', '.fasta'), 'wt') as out:
        for i, s in output:
            out.write('>' + i + '\n' + s + '\n')
"

combine_stats="#!/usr/bin/env python3
import csv
import itertools
import os
import re
import sys
from Bio import SeqIO
import numpy

replicons = {}
with open('${REPLICONS}', 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        replicons[row[0]] = row[6]

gwa_strains = set()
with open('${GWA_STRAINS}', 'rt') as handle:
    for line in handle: gwa_strains.add(line.strip())

data = {}
all_stats = set()
with open('compute_output.tsv', 'rt') as handle:
    compute_stats = []
    for line in handle:
        line = line.strip()
        if line.startswith('#') or len(line) == 0: continue
        if line.startswith('locus'):
            compute_stats = line.split('\t')
            continue
        else:
            stats = line.split('\t')
            gene = re.sub('\.fa[a-z]+', '', stats[0]).split('/')[1]
            data[gene] = {}
            for cs, s in zip(compute_stats[1:], stats[1:]):
                data[gene][cs] = s
        seqfname = 'protein/' + gene + '.fasta'
        if not os.access(seqfname, os.F_OK):
            seqfname = 'nucleotide/' + gene + '.fasta'
        strains = set(); n_seqs = 0
        g_strains = set(); g_n_seqs = 0
        for rec in SeqIO.parse(seqfname, 'fasta'):
            strain = rec.id.split('__')[0]
            strains.add(strain)
            n_seqs += 1
            if strain in gwa_strains:
                g_strains.add(strain)
                g_n_seqs += 1
        data[gene]['n_strains'] = str(len(strains))
        data[gene]['n_sequences'] = str(n_seqs)
        data[gene]['gwa_n_strains'] = str(len(g_strains))
        data[gene]['gwa_n_sequences'] = str(g_n_seqs)
        all_stats.update(data[gene].keys())

all_stats.update(['beta_mean', 'beta_median', 'beta_max'])
with open('${BETA}', 'rt') as handle:
    for gene, rows in itertools.groupby(csv.reader(handle, delimiter=' '), lambda r: r[0]):
        values = [float(r[2]) for r in rows]
        data[gene]['beta_mean'] = numpy.mean(values)
        data[gene]['beta_median'] = numpy.median(values)
        data[gene]['beta_max'] = max(values)

all_stats.update(['raisd_mean', 'raisd_median', 'raisd_max'])
with open('${RAISD}', 'rt') as handle:
    for gene, rows in itertools.groupby(csv.reader(handle, delimiter=' '), lambda r: r[0]):
        values = [float(r[2]) for r in rows]
        data[gene]['raisd_mean'] = numpy.mean(values)
        data[gene]['raisd_median'] = numpy.median(values)
        data[gene]['raisd_max'] = max(values)


with open('pdnds_output.tsv', 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        gene = re.sub('\.fa[a-z]+', '', row[0]).split('/')[1]
        i = 0
        for sites in ['all', 'exons', 'introns', 'nonsyn', 'syn', 'thirds', 'fourfold', 'silent']:
            for stat in ['S', 'n_mutations', 'n_singletons', 'n_sites', 'ThetaW', 'ThetaPi', 'TajD']:
                i += 1
                if sites in {'nonsyn', 'syn', 'fourfold'}:
                    data[gene][sites + '_' + stat] = row[i]
        all_stats.update(data[gene].keys())
 
with open('maf.tsv', 'rt') as handle:
    all_stats.update(['mean_maf', 'median_maf'])
    for row in csv.DictReader(handle, delimiter='\t'):
        if row['maf'] != '':
            maf = list(map(float, row['maf'].split(',')))
            g = row['gene'].replace('.vcf', '')
            if len(maf) > 0:
                data[g]['mean_maf'] = numpy.mean(maf)
                data[g]['median_maf'] = numpy.median(maf)
            else:
                data[g]['mean_maf'] = float('nan')
                data[g]['median_maf'] = float('nan')

output_stats = ['gene', 'replicon'] + list(all_stats - {'gene'})
print('\t'.join(output_stats))
for gene, row in data.items():
    sys.stdout.write(gene)
    for s in output_stats:
        if s == 'gene': continue
        if s == 'replicon':
            sys.stdout.write('\t' + replicons[gene])
        elif s in row:
            sys.stdout.write('\t' + str(row[s]))
        else:
            sys.stdout.write('\tnan')
    sys.stdout.write('\n')
"


if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load analysis/0.8.8
    module load gcc
    module load R/4.0.4

    set -euo pipefail

    SECONDS=0

    cd "$OUTDIR"

    TASK="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    PHENOTYPE="$(cut -f 2 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    case "$TASK" in

        prep)
            echo "$trim_and_filter" > trim_and_filter.py; chmod u+x trim_and_filter.py
            awk -F$'\t' 'NR == 1 { print };' "$ORTHOLOGS" \
                | tr '\t' '\n' \
                | grep magncbi \
                | sed 's/__.\+//g' \
                | sort -u \
                > strains.txt \
                || { echo "making list of strains failed"; exit 1; }

            ## Copy alignments while trimming off terminal Ns and gaps and
            ## filter to just the strains of interest
            rm -rf protein nucleotide
            mkdir -p protein nucleotide
            ./trim_and_filter.py "$ALIGNMENT_DIRECTORY/aligned-protein" strains.txt protein \
                || { echo "filtering protein alignments failed"; exit 1; }
            ./trim_and_filter.py "$ALIGNMENT_DIRECTORY/aligned-nucleotide" strains.txt nucleotide \
                || { echo "filtering nucleotide alignments failed"; exit 1; }
            ;;

        calc-popgen)
            ## Run compute (theta, etc.) and polydNdS (ns and syn) on
            ## every gene
            ## Does not remove ambiguous sites, but does remove any site
            ## with a deletion.
            rm -f commands
            for fname in protein/*.fasta nucleotide/*.fasta; do
                echo "compute -i \"${fname}\"" >> commands
            done
            singularity shell $ANALYSIS_DIR/lsanalysis < commands \
                > compute_output.tsv \
                || { echo "running compute failed"; exit 1; }

            rm -f commands
            for fname in protein/*.fasta; do
                # -T: parseable output
                # -A: approximate treatment of highly diverged codons
                echo "polydNdS -i \"${fname}\" -T -A" >> commands
            done
            singularity shell $ANALYSIS_DIR/lsanalysis < commands \
                > pdnds_output.tsv \
                || { echo "running polydNdS failed"; exit 1; }

            ## Minor allele frequency
            python3 -c "$maf" "$VCF" > maf.tsv \
                || { echo "calculating maf failed"; exit 1; }

            ## gestimator kaks
            ## Added 2021-08-02
            rm -f commands
            for fname in "${BBH_ALIGNS}protein/"*.fas; do
                gene="$(basename "${fname/.fas/}")"
                echo "gestimator -i \"${fname}\" | awk '{print \"${gene}\t\" "'$0'"};'" >> commands
            done
            echo "medmel_ortholog	gene1	gene2	ka	ks	ka_ks" > "gestimator_output.tsv"
            singularity shell $ANALYSIS_DIR/lsanalysis < commands \
                | sed 's/\t999\t/\tNaN\t/g' \
                | sed 's/\t999$/\tNaN/g' \
                >> gestimator_output.tsv \
                || { echo "calculating KaKs failed"; exit 1; }

            ## Combine with a table of counts: strains, sequences, copy number
            echo "$combine_stats" > combiner.py; chmod u+x combiner.py
            ./combiner.py > stats.tsv \
                || { echo "combining tables failed"; exit 1; }
            ;;

        nssites)
            ## Get a list of synonymous and non-synonymous sites for every gene
            rm -f nsyn_sites.txt syn_sites.txt 4f_sites.txt
            for fname in protein/*.fasta; do
                gene="$(basename "${fname/.fasta/}")"
                ## Change gaps to N so gapped sites won't be skipped
                sed 's/-/N/g' "$fname" > tmp \
                    || { echo "changing gaps to N for ${fname} failed"; exit 1; }
                singularity shell "${ANALYSIS_DIR}/lsanalysis" \
                    <<< "polydNdS -P -A -i tmp > /dev/null" \
                    || { echo "running on ${fname} failed"; exit 1; }
                awk -F$'\t' 'NR == 2 { for(i=1; i < NF; i++) { print "'$gene'\t" $i }; };' \
                    tmp.replacement \
                    >> nsyn_sites.txt \
                    || { echo "print non-synonymous sites for ${fname} failed"; exit 1; }
                awk -F$'\t' 'NR == 2 { for(i=1; i < NF; i++) { print "'$gene'\t" $i }; };' \
                    tmp.synonymous \
                    >> syn_sites.txt \
                    || { echo "print synonymous sites for ${fname} failed"; exit 1; }
                awk -F$'\t' 'NR == 2 { for(i=1; i < NF; i++) { print "'$gene'\t" $i }; };' \
                    tmp.4fold \
                    >> 4f_sites.txt \
                    || { echo "print 4fold sites for ${fname} failed"; exit 1; }
            done
            rm tmp

            # Table of MAF by variant, split into biallelic sites
            ## Use VCF2 to make sure alternate states are coded the same
            ## way as the association analyses
            echo "$maf_table" > maf_table.py && chmod u+x maf_table.py
            ./maf_table.py "$VCF2" \
                > variant_af_table.tsv
            ;;


    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"

    rm "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-${1}-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    case "$1" in
        prep)
            WALLTIME=00:30:00
            NTHREADS=1
            MEM=2GB
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        calc-popgen)
            WALLTIME=16:00:00
            NTHREADS=1
            MEM=64GB
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        nssites)
            WALLTIME=10:00:00
            NTHREADS=1
            MEM=8GB
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --ntasks-per-node=$NTHREADS \
        --mem=$MEM --time=$WALLTIME --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
