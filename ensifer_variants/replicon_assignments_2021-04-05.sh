#!/bin/bash
#
# Assign replicons to orthologs from the 2020-04-09 whole genome alignment.
#
# Updated with new orthologs 2021-04-05.
#
# Updated 2021-06-15 with "filtered" dataset and code to calculate
# the median position in the replicon.
#
# SETTINGS
#
QUEUE=amdsmall
WALLTIME="00:05:00"
MEM=8GB
NTHREADS=1
#

fai_to_replicon="#!/usr/bin/env python3
import csv
import sys
for fname in sys.argv[1:]:
    lengths = []
    strain = fname.split('.', 1)[0]
    with open(fname, 'rt') as handle:
        for row in csv.reader(handle, delimiter='\t'):
            lengths.append((int(row[1]), row[0]))
    if len(lengths) <= 7:
        lengths.sort(reverse=True)
        for i, (l, c) in enumerate(lengths):
            if i == 0: print(c, 'chromosome', l, sep='\t')
            elif i == 1: print(c, 'psymb', l, sep='\t')
            elif i == 2: print(c, 'psyma', l, sep='\t')
            else: print(c, 'plasmid', l, sep='\t')
"

assign_replicons="#!/usr/bin/env python3
import collections
import csv
import gzip
import re
import sys

import numpy

def openfun(fname, mode):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

contigs = {}
lengths = {}
with openfun(sys.argv[1], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        contigs[re.sub('__.+?\.', '.', row[0])] = row[1]
        lengths[re.sub('__.+?\.', '.', row[0])] = int(row[2])

## Because the annotation table does not include assembly versions,
## there is a small chance of confusion among replicons. In practice,
## I don't think this will occur.
genes = collections.defaultdict(lambda: {'chromosome':0, 'psyma':0, 'psymb':0, 'plasmid':0, 'total':0})
positions = collections.defaultdict(list)
with openfun(sys.argv[2], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        ol, strain, _, contig, lt = row[:5]
        s, e = tuple(map(int, row[8:10]))
        ctg = strain + '.' + contig
        if ctg in contigs:
            genes[ol]['total'] += 1
            genes[ol][contigs[ctg]] += 1
            positions[ol].append(((s+e)/2)/lengths[ctg])
        else:
            _ = genes[ol]['total']

rep_order = ('chromosome', 'psymb', 'psyma', 'plasmid')
for ol in sorted(genes, key=lambda x: int(x.replace('OL', ''))):
    counts = genes[ol]
    total = genes[ol]['total']
    median_pos = numpy.median(positions[ol])
    if total > 0:
        rep_counts = (counts[rep_order[0]]/total, counts[rep_order[1]]/total,
                      counts[rep_order[2]]/total, counts[rep_order[3]]/total)
        if max(rep_counts) > 0.5:
            assignment = rep_order[rep_counts.index(max(rep_counts))]
            pl_assignment = assignment
            if pl_assignment == 'plasmid': pl_assignment = 'psyma'
        else:
            assignment = 'unclassified'
            if rep_counts[rep_order.index('psyma')] + rep_counts[rep_order.index('plasmid')] > 0.5:
                pl_assignment = 'psyma'
            else:
                pl_assignment = 'unclassified'
    else:
        rep_counts = (0, 0, 0, 0)
        assignment = 'unclassified'
        pl_assignment = 'unclassified'
    print(ol, '\t'.join(map(str, rep_counts)), assignment, pl_assignment, median_pos, sep='\t')
"

PROJDIR="/home/tiffinp/epste051/project/jgi_sequencing"
INDIR="${PROJDIR}/results/ortho/medmel/sibeliaz/2020-04-09_split"

OUTDIR="${PROJDIR}/results/replicon_assignments/2021-04-05"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    set -euo pipefail

    SECONDS=0

    SPECIES="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    cd "$OUTDIR"
    mkdir -p "$SPECIES" && cd "$SPECIES"

    echo "$fai_to_replicon" > fai_to_replicon.py; chmod u+x fai_to_replicon.py
    ./fai_to_replicon.py "${INDIR}/${SPECIES}/fasta/"*.fai \
        > contigs_to_replicons.tsv \
        || { echo "getting replicon names for ${SPECIES} failed"; exit 1; }

    echo "$assign_replicons" > assign_replicons.py; chmod u+x assign_replicons.py
    ./assign_replicons.py contigs_to_replicons.tsv \
        "${INDIR}/${SPECIES}/annotation_table.overlap.0.80.filtered.tsv.gz" \
        > replicon_assignments.filtered.tsv \
        || { echo "assigning replicons failed"; exit 1; }
        # Previous input file: #"${INDIR}/${SPECIES}/annotation_table.overlap.0.80.tsv.gz" \
        #> replicon_assignments.filtered.tsv \

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
    for species in medicae meliloti; do
        echo "$species" > "${ARRAYDIR}/${i}"
        i=$(($i+1))
    done

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition="$QUEUE" --time=$WALLTIME --mem=$MEM \
        --nodes=1 --ntasks=$NTHREADS --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
