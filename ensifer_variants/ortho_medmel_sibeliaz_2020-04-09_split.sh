#!/bin/bash
#
# Run SibeliaZ to get whole genome alignments of all
# the medicae and meliloti strains separately.
#
# OUTPUT FILES (per species)
#
# Note that some of these files may be gzipped or deleted if they are very
# large.
#
#   all.maf.gz                      Whole genome alignment of all LCBs greater than the minimum size.
#   annotation.tsv.gz               Features overlapping alignment blocks in output/all_blocks.gff. Due to contig naming issues, will not have all versions of the annotations.
#   annot.gff3                      Single concatenated file with annotations for all strains. (Might be deleted or gzipped.)
#   jaccard.size_filtered.tsv       Jaccard index for all pairs of strains based on sharing of LCBs that have been filtered by minimum size.
#   lcbs.size_filtered.vcf.gz       VCF file of LCB presence/absence (not copy number). Only LCBs > a minimum size.
#   lcb_table.size_filtered.tsv     Tsv file with LCB copy number for each strain + some length information. Only LCBs > a minimum size.
#   lcb_table.tsv                   From custom script. Tsv file with LCB copy number for all LCBs (in output/all_blocks.gff).
#   output/all_blocks.gff           From custom script. Combination of aligned and unaligned regions (see next two files)
#   output/blocks_coords.gff        From SibeliaZ. Contains coordinates of aligned regions.
#   output/unaligned.gff            From custom script. Contains coordinates of regions that didn't align to other strains.
#   resource_usage.txt              Wall time and memory usage for the SibeliaZ alignment search step.
#
# SETTINGS
#
QUEUE=amdsmall
WALLTIME="06:00:00"
MEM=64GB
NTHREADS=64
OLDPYTHONPATH="$HOME/usr/lib/python2.7/site-packages.old"
#
MAX_BUBBLE_SIZE=1000    # -b parameter, distance between points in a graph bubble; default=200
MAX_COPIES_SCALAR=10    # -a parameter after multiplying by the number of genomes, maximum number of copies of kmer allowed per genome
MIN_BLOCK_SIZE=50       # -m parameter; minimum LCB size; default=50
KMER=15                 # -k parameter; kmer size for searches; default=25, 15 recommended for small datasets, makes little difference for close datasets (according to authors)
N_CHUNKS=200            # number of alignment processes to run
#
N_OVERLAP_CHUNKS=300

declare -A STRAIN_PREFIXES
STRAIN_PREFIXES["complete/others"]="crfs"
STRAIN_PREFIXES["complete/mage_pacbio_2016"]="pbmg"
STRAIN_PREFIXES["draft/mage_ensifer_2010"]="inr"
STRAIN_PREFIXES["draft/new_short_read_2016"]="inb"
STRAIN_PREFIXES["draft/mag_strains_ncbi"]="magncbi"

declare -A EXTRA_STRAIN_PREFIXES
EXTRA_STRAIN_PREFIXES["draft/mag_strains_rast"]="magrast"
EXTRA_STRAIN_PREFIXES["complete/ncbi_pacbio_2016"]="pbnb"
#
BAD_STRAINS=("5A14__drfs" "A0641M__drfs" "JGI0001011_A08__drfs" "PC2__drfs" "RP12G__drfs") # Known poor assemblies to avoid
#
SPECIES_LIST=(medicae meliloti)

gff_filter="
import os
import sys
import tempfile

import gffutils

tf = tempfile.mkstemp(dir='.')[1]
genes = {}
n = 0
with open(tf, 'wt') as out:
    done = set()
    with open(sys.argv[1], 'rt') as handle:
        for line in handle:
            if line.startswith('#'): continue
            row = line.strip().split('\t')
            coords = (row[3], row[4])
            try:
                ID = [x.split('=',1)[1] for x in row[8].split(';') if x.startswith('ID=')][0]
            except IndexError:
                ## Line does not have an ID; some old MaGe file
                if row[2] == 'gene':
                    ID = 'gene' + str(n)
                    n += 1
                    if ID in done:
                        raise Exception('Can\'t get a unique ID: %s' % row[8])
                    row[8] = 'ID=' + ID + ';' + row[8]
                    genes[coords] = ID
                else:
                    ## Not a gene and lacking an ID, keep going
                    continue
            if row[2] != 'gene' and 'Parent=' not in row[8]:
                if coords in genes:
                    row[8] = row[8] + ';Parent=' + genes[coords]
            if ID not in done:
                out.write('\t'.join(row) + '\n')
                done.add(ID)

def fill_info(feature, data):
    if 'locus_tag' in feature.attributes:
        if isinstance(feature.attributes['locus_tag'], list):
            ## If there is more than one, we have a problem, but I
            ## don't think that happens in our dataset
            lt = feature.attributes['locus_tag'][0]
        else:
            lt = feature.attributes['locus_tag']
        if data['locus_tag'] is not None and lt != data['locus_tag']:
            raise Exception('locus_tag does not match %s, %s' % (lt, data['locus_tag']))
        data['locus_tag'] = lt
    if 'pseudo' in feature.attributes and not data['pseudo']:
        data['pseudo'] = True
    if feature.featuretype == 'pseudogene' and not data['pseudo']:
        data['pseudo'] = True
    if 'partial' in feature.attributes and not data['fragment']:
        data['fragment'] = True
    if 'gene_biotype' in feature.attributes and feature.attributes['gene_biotype'] == 'pseudogene':
        data['pseudo'] = True
    if 'product' in feature.attributes:
        for p in feature.attributes['product']:
            data['product'].add(p)
            if 'fragment' in p:
                data['fragment'] = True

db = gffutils.create_db(tf, ':memory:', merge_strategy='merge')
for flist in db.iter_by_parent_childs(featuretype='gene'):
    data = {'locus_tag':None, 'pseudo':False, 'fragment':False, 'product':set()}
    parent = flist[0]
    fill_info(parent, data)
    if len(flist) > 1:
        for child in flist[1:]:
            fill_info(child, data)
            for child2 in db.children(child.id):
                fill_info(child2, data)
    if data['fragment'] or data['pseudo']: continue
    output = [parent.chrom, parent.source, parent.featuretype,
              parent.start, parent.end, parent.score, parent.strand,
              '0',
              'ID=' + parent.attributes['ID'][0] + ';' +
              'locus_tag=' + data['locus_tag'] + ';' +
              'product=' + ','.join(data['product'])]
    sys.stdout.write('\t'.join(map(str, output)) + '\n')

## Need to keep the regions for some downstream steps
for flist in db.iter_by_parent_childs(featuretype='region'):
    parent = flist[0]
    output = [parent.chrom, parent.source, parent.featuretype,
              parent.start, parent.end, parent.score, parent.strand,
              '0', '.']
    sys.stdout.write('\t'.join(map(str, output)) + '\n')

os.unlink(tf)
"

change_replicon_names="
import csv
import sys
to_change = []
rdr = csv.reader(sys.stdin, delimiter='\t')
to_change_reps = {}
for row in rdr:
    if not row[0].startswith('#'):
        if row[2] == 'region':
            to_change_reps[row[0]] = row[4]
    to_change.append(row)
target_reps = {}
with open(sys.argv[1], 'rt') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        if not row[0].startswith('#'):
            if row[2] == 'remark':
                if row[4] in list(to_change_reps.values()):
                    target_reps[row[4]] = row[0]
                else:
                    target_reps[str(int(row[4])-1)] = row[0]
sys.stderr.write(str(target_reps) + '\n')
sys.stderr.write(str(to_change_reps) + '\n' + '\n')
for row in to_change:
    if not row[0].startswith('#'):
        row[0] = target_reps[to_change_reps[row[0]]]
    print('\t'.join(row))
"

extract_locus_tag="
import csv
import re
import sys
rdr = csv.reader(sys.stdin, delimiter='\t')
for row in rdr:
    a = {x.split('=',1)[0]: x.split('=',1)[1] for x in re.sub(';$', '', row[3]).split(';')}
    if 'locus_tag' in a:
        row[3] = a['locus_tag']
    else:
        continue
    #elif 'ID' in a:
    #    row[3] = a['ID']
    print('\t'.join(row))
"

annotation_table="
import csv
import gzip
import sys
with gzip.open(sys.argv[1], 'rt') as handle:
    locus_tags = {}
    for line in handle:
        row = line.strip().split('\t')
        strain, contig = row[0].split('.', 1)
        strain = strain.split('__')[0]
        lcb = row[8].replace('id=', '')
        asource = ''
        lt = ''; p = ''; ps = ''
        if len(row) >= 18 and row[17] != '.':
            if row[10] == 'RAST': asource = 'rast'
            if row[10] == 'MAGE': asource = 'mage'
            if row[10] == 'NCBI': asource = 'ncbi'
            tags = {x.split('=',1)[0]: x.split('=',1)[1] for x in row[17].split(';')}
            ps = '0'
            if 'locus_tag' in tags:
                lt = tags['locus_tag']
                locus_tags[(row[12], row[13])] = lt
            elif 'Name' in tags and ':' in tags['Name']:
                lt = tags['Name'].split(':', 1)[1]
            elif 'ID' in tags:
                lt = tags['ID']
            if 'product' in tags: p = tags['product']
            if 'locus_tag' not in tags and 'product' in tags:
                if (row[12], row[13]) in locus_tags:
                    lt = locus_tags[(row[12], row[13])]
            if 'pseudo' in tags: ps = '1'
            if p == '': continue
        print(lcb, strain, asource, contig, lt, p, ps, sep='\t')
"

annotation_table_from_ologs="#!/usr/bin/env python3
import csv
import gzip
import itertools
import sys
ologs = {}
with open(sys.argv[1], 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        for k, v in row.items():
            if k == 'ortholog':
                o = row[k]
            else:
                for g in row[k].strip().split(','):
                    if g == '': continue
                    ologs[(k.split('__')[0], g)] = o
with gzip.open(sys.argv[2], 'rt') as handle:
    done = set()
    parent_gene_tuples = {}
    locus_tags = {}
    partials = set()
    pseudos = set()
    products = {}
    for gene_tuple, rows in itertools.groupby(csv.reader(handle, delimiter='\t'), lambda r: (r[9], r[12], r[13])):
        lt = ''; p = ''; ps = '0'; pt = '0'
        for row in rows:
            strain, contig = row[0].split('.', 1)
            strain = strain.split('__')[0]
            asource = ''
            start, end = row[12], row[13]
            if len(row) >= 18 and row[17] != '.':
                if row[10] == 'RAST': asource = 'rast'
                if row[10] == 'MAGE': asource = 'mage'
                if row[10] == 'NCBI': asource = 'ncbi'
                tags = {x.split('=',1)[0]: x.split('=',1)[1] for x in row[17].split(';')}
                if 'Parent' in tags and (row[0], tags['Parent']) in parent_gene_tuples:
                    gene_tuple = parent_gene_tuples[(row[0], tags['Parent'])]
                if 'ID' in tags:
                    parent_gene_tuples[(row[0], tags['ID'])] = gene_tuple
                if 'locus_tag' in tags:
                    lt = tags['locus_tag']
                    locus_tags[gene_tuple] = lt
                if 'product' in tags:
                    p = tags['product']
                    products[gene_tuple] = p
                if 'locus_tag' not in tags and 'product' in tags:
                    if gene_tuple in locus_tags:
                        lt = locus_tags[gene_tuple]
                if lt != '' and p == '' and gene_tuple in products:
                    p = products[gene_tuple]
                if 'pseudo' in tags or gene_tuple in pseudos:
                    ps = '1'
                    pseudos.add(gene_tuple)
                if 'partial' in tags or gene_tuple in partials:
                    partials.add(gene_tuple)
                    pt = '1'
                if 'product' in tags and 'fragment' in tags['product']:
                    partials.add(gene_tuple)
                    pt = '1'
                if p == '': continue
            if lt == '' or p == '': continue
            try:
                o = ologs[(strain, lt)]
            except KeyError:
                continue
        if (o, strain, lt) in done:
            continue
        if (strain, lt) not in ologs:
            continue
        if p == '':
            continue
        done.add((o, strain, lt))
        print(o, strain, asource, contig, lt, p, ps, pt, start, end, sep='\t')
"

print_orthologs="
import collections
import os
import sys
import igraph
with open(sys.argv[2], 'rt') as handle:
    strains = sorted(l.strip() for l in handle.readlines())
cc = igraph.load(sys.argv[1], format='ncol', directed=False).clusters(mode=igraph.WEAK)
sys.stdout.write('ortholog\t' + '\t'.join(strains) + '\n')
all_genes = set()
for o, c in enumerate(cc):
    sys.stdout.write('OL' + str(o))
    genes = collections.defaultdict(set)
    for n in c:
        s, g = cc.graph.vs[n]['name'].split(':')
        genes[s].add(g)
        all_genes.add(s + ':' + g)
    for strain in strains:
        if strain in genes:
            sys.stdout.write('\t' + ','.join(genes[strain]))
        else:
            sys.stdout.write('\t')
    sys.stdout.write('\n')
## This block commented out 2021-06-03:
'''
for fname in os.listdir(sys.argv[3]):
    strain = fname.split('.')[0]
    if not strain in strains: continue
    with open(sys.argv[3] + '/' + fname, 'rt') as handle:
        for line in handle:
            row = line.strip().split('\t')
            if row[2] != 'gene': continue
            tags = {x.split('=',1)[0]: x.split('=',1)[1] for x in row[8].split(';')}
            if 'gene_biotype' in tags and tags['gene_biotype'] == 'pseudogene': continue
            if 'pseudo' in tags: continue
            if 'locus_tag' in tags:
                lt = tags['locus_tag']
                if strain + ':' + lt not in all_genes:
                    o += 1
                    sys.stdout.write('OL' + str(o))
                    for s in strains:
                        if s == strain:
                            sys.stdout.write('\t' + lt)
                        else:
                            sys.stdout.write('\t')
                    sys.stdout.write('\n')
'''
"

PROJDIR="/home/tiffinp/epste051/project/jgi_sequencing"
ASSMDIR="${PROJDIR}/data/assembly"
KEEP="${PROJDIR}/data/metadata/strain_versions_to_keep.2019-12-27.txt"

OUTDIR="${PROJDIR}/results/ortho/medmel/sibeliaz/2020-04-09_split"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/ortho_2020-04-09_split"

unaligned="${PROJDIR}/script/bin/process_sibeliaz_gff.py"
tabler="${PROJDIR}/script/bin/annotate_sibeliaz_blocks.py"
vcf="${PROJDIR}/script/bin/lcbtsv2vcf.py"
jacc="${PROJDIR}/script/bin/jaccard_index.py"
maf_aligner="${PROJDIR}/script/bin/run_spoa_maf.py"
varfinder="${PROJDIR}/script/bin/variants_from_maf.py"
make_order_files="${PROJDIR}/script/bin/sibeliaz_gff_to_lcb_order.py"
compare_order="${PROJDIR}/script/bin/comp_order"
coord_table="${PROJDIR}/script/bin/maf_to_coords_table.py"
lcb_len="${PROJDIR}/script/bin/maf_lcb_len.py"
tally_overlap="${PROJDIR}/script/bin/overlap_tally_by_strain.py"
overlap_tally_to_ologs="${PROJDIR}/script/bin/overlap_tally_to_orthologs2.py"
make_links="${PROJDIR}/script/bin/sorted_coords_to_links.py"
grab_sequence="${PROJDIR}/one_sequence_per_lcb_gfa.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools
    module load bedtools/2.29.0
    module load htslib/1.9
    module load mummer/4.0.0.beta2
    module load parallel/20190122
    module load parsnp/1.5.6
    module load R/3.4.3
    module load sibeliaz/1.2.0
    module load spoa/3.0.1
    module load vg/1.29.0

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    SPECIES="$(cut -f 2 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    SUBDIR="$(cut -f 3 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    CHUNK="$(cut -f 4 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    cd "$OUTDIR"
    mkdir -p "$SPECIES" && cd "$SPECIES"

    case "$TASK" in

        "get-sequences")

            mkdir -p "fasta"
            mkdir -p "gff"

            group="${STRAIN_PREFIXES[$SUBDIR]}"
            for f in "${ASSMDIR}/${SUBDIR}/"*; do
                strain="$(basename "$f")" ## Will not work with certain whitespace symbols in filenames
                if [[ "$strain" =~ "README" || "$strain" =~ "." ]]; then
                    continue
                fi

                ## Get rid of strains that are not listed as the correct species
                grep "Sinorhizobium ${SPECIES}" \
                    "${ASSMDIR}/${SUBDIR}/${strain}/species.txt" || \
                    continue

                ## We don't want any dashes in the names
                new_strain_id="$(echo "$strain" | sed 's/-/_/g')"
                echo "Working on ${strain}__${group} at ${SECONDS}"


                echo "    writing ${ASSMDIR}/${SUBDIR}/${strain}/genome.fasta to fasta/${new_strain_id}__${group}.fasta"
                sed 's/^>/>'"${new_strain_id}__${group}"'./g' \
                    "${ASSMDIR}/${SUBDIR}/${strain}/genome.fasta" > \
                    "fasta/${new_strain_id}__${group}.fasta" \
                    || { echo "copying fasta ${new_strain_id}__${group} failed"; exit 1; }

                echo "    writing ${ASSMDIR}/${SUBDIR}/${strain}/genome.gff3 to gff/${new_strain_id}__${group}.gff3"
                grep -v "^#" "${ASSMDIR}/${SUBDIR}/${strain}/genome.gff3" \
                    | sed 's/^/'"${new_strain_id}__${group}"'./g' \
                    | sed 's/\.NODE/.lcl|NODE/g' \
                    > "gff/${new_strain_id}__${group}.gff3" \
                    ||  { echo "copying gff ${new_strain_id}__${group} failed"; exit 1; }

            done

            ## Remove the "bad" strains
            for s in "${BAD_STRAINS[@]}"; do
                rm -f "fasta/${s}."* \
                    || { echo "deleting fasta for ${s} failed"; exit 1; }
                rm -f "gff/${s}."* \
                    || { echo "deleting gff for ${s} failed"; exit 1; }
            done
            ;;

        "align")

            n_strains="$(ls "fasta" | grep -c "fasta$")"
            max_copy=$(($MAX_COPIES_SCALAR*$n_strains))

            /usr/bin/time -v -o "resource_usage.${KMER}.${MAX_BUBBLE_SIZE}.txt" \
                sibeliaz -n -o "output.${KMER}.${MAX_BUBBLE_SIZE}" \
                -t "$NTHREADS" \
                -k "$KMER" \
                -m "$MIN_BLOCK_SIZE" \
                -a "$max_copy" \
                -b "$MAX_BUBBLE_SIZE" \
                fasta/*.fasta \
                || { echo "aligning failed"; exit 1; }
            ;;

        "process")

#            ## Add in different annotation versions
#            ## UPDATE 2020-12-21: Fixed replicon names so this extra
#            ## annotation information will actually work.
#            subdir="complete/ncbi_pacbio_2016"
#            group="${EXTRA_STRAIN_PREFIXES[$subdir]}"
#            for f in "${ASSMDIR}/${subdir}/"*; do
#                strain="$(basename "$f")" ## Will not work with certain whitespace symbols in filenames
#                if [[ "$strain" =~ "README" || "$strain" =~ "." ]]; then
#                    continue
#                fi
#
#                ## Get rid of strains that are not listed as medicae or meliloti
#                grep "Sinorhizobium ${SPECIES}" \
#                    "${ASSMDIR}/${subdir}/${strain}/species.txt" || \
#                    continue
#
#                ## We don't want any dashes in the names
#                new_strain_id="$(echo "$strain" | sed 's/-/_/g')"
#
#                ## Modify file to match other annotations
#                grep -v "^#" "${ASSMDIR}/${subdir}/${strain}/genome.gff3" \
#                    | sed 's/^/'"$new_strain_id"'__pbmg./g' \
#                    | python -c "$change_replicon_names" "gff/${new_strain_id}__pbmg.gff3" \
#                    | awk -F$'\t' '{OFS="\t"; $2="NCBI"; print}' \
#                    > "gff/${new_strain_id}__${group}.gff3" \
#                    ||  { echo "copying gff ${new_strain_id}__${group} failed"; exit 1; }
#            done
#
#            subdir="draft/mag_strains_rast"
#            group="${EXTRA_STRAIN_PREFIXES[$subdir]}"
#            for f in "${ASSMDIR}/${subdir}/"*; do
#                strain="$(basename "$f")" ## Will not work with certain whitespace symbols in filenames
#                if [[ "$strain" =~ "README" || "$strain" =~ "." ]]; then
#                    continue
#                fi
#
#                ## Get rid of strains that are not listed as medicae or meliloti
#                grep "Sinorhizobium meliloti\|Sinorhizobium medicae" \
#                    "${ASSMDIR}/${subdir}/${strain}/species.txt" || \
#                    continue
#
#                ## We don't want any dashes in the names
#                new_strain_id="$(echo "$strain" | sed 's/-/_/g')"
#
#                ## Modify file to match other annotations
#                grep -v "^#" "${ASSMDIR}/${subdir}/${strain}/genome.gff3" \
#                    | sed 's/^/'"$new_strain_id"'__magncbi./g' \
#                    | sed 's/\.NODE/.lcl|NODE/g' \
#                    | awk -F$'\t' '{OFS="\t"; $2="RAST"; print}' \
#                    > "gff/${new_strain_id}__${group}.gff3" \
#                    ||  { echo "copying gff ${new_strain_id}__${group} failed"; exit 1; }
#            done
#
#
#            ## Add in unaligned sequences
#            cp "$unaligned" .
#            "./$(basename "$unaligned")" --output "output.${KMER}.${MAX_BUBBLE_SIZE}/unaligned.gff" \
#                "output.${KMER}.${MAX_BUBBLE_SIZE}/blocks_coords.gff" "fasta" \
#                || { echo "finding unaligned blocks failed"; exit 1; }
#            cat "output.${KMER}.${MAX_BUBBLE_SIZE}/blocks_coords.gff" "output.${KMER}.${MAX_BUBBLE_SIZE}/unaligned.gff" > \
#                "output.${KMER}.${MAX_BUBBLE_SIZE}/all_blocks.gff" \
#                || { echo "combining blocks failed"; exit 1; }

            ## Annotations.  Intersect with output.
            cat gff/*.gff3 \
                > "annot.gff3" \
                || { echo "combining annotations failed"; exit 1; }
            bedtools intersect -a "output.${KMER}.${MAX_BUBBLE_SIZE}/all_blocks.gff" \
                -b "annot.gff3" \
                -wao \
                > "annotation.tsv" \
                || { echo "annotating ${SPECIES} failed"; exit 1; }
            gzip -f annotation.tsv

            ## ADDED 2020-12-21
            ## Make a more readable table
            python -c "$annotation_table" "annotation.tsv.gz" \
                > annotation_table.tsv \
                || { echo "making annotation table for ${SPECIES} failed"; exit 1; }
            gzip -f annotation_table.tsv

#            ## Table of LCB copies
#            cp "$tabler" .
#            "./$(basename "$tabler")" "output.${KMER}.${MAX_BUBBLE_SIZE}/all_blocks.gff" > \
#                "lcb_table.tsv" \
#                || { echo "making table of LCB copy number failed"; exit 1; }
#            awk -F$'\t' '$5 >= '"$MIN_BLOCK_SIZE"' {print};' "lcb_table.tsv" > \
#                "lcb_table.size_filtered.tsv" \
#                || { echo "filtering on LCB size failed"; exit 1; }
#
#            cp "$vcf" .
#            "./$(basename "$vcf")" --keep "$KEEP" "lcb_table.size_filtered.tsv" | \
#                sed 's/__.\{3,7\}\t/\t/g' | sed 's/__.\{3,7\}$//g' > \
#                "lcbs.size_filtered.vcf" \
#                || { echo "converting to VCF failed"; exit 1; }
#            rm -f "lcbs.size_filtered.vcf.gz"
#            bgzip -i "lcbs.size_filtered.vcf" \
#                || { echo "bgzipping failed"; exit 1; }
#            tabix "lcbs.size_filtered.vcf.gz" \
#                || { echo "tabixxing failed"; exit 1; }
#
#            cp "$jacc" .
#            "./$(basename "$jacc")" --vcf "lcbs.size_filtered.vcf.gz" "_" > \
#                "jaccard.size_filtered.tsv" \
#                || { echo "calculating jaccard index failed"; exit 1; }
            ;;

        "maf")
            # #cat fasta/*.fasta > all.fasta ## Do manually
            # samtools faidx all.fasta

            echo "listing LCBs for chunk ${CHUNK}"
            cut -f 1 "lcb_table.size_filtered.tsv" | \
                tail -n +2 | \
                split -n "r/${CHUNK}/${N_CHUNKS}" > \
                "lcbs_to_align.${CHUNK}.txt" \
                || { echo "getting list of LCBs to align failed"; exit 1; }

            mkdir -p "alignments"
            rm -f "alignments/${CHUNK}.maf"
            mkdir -p "gff_chunks_${CHUNK}"
            mkdir -p "fasta_chunks_${CHUNK}"
            while read -r lcb; do
                ## Extract fasta sequence...
                ## bedtools getfasta from concatenated sequences
                echo "LCB=${lcb}, chunk=${CHUNK}"
                grep "id=${lcb}$" "output.${KMER}.${MAX_BUBBLE_SIZE}/all_blocks.gff" > \
                    "gff_chunks_${CHUNK}/${lcb}.gff" \
                    || { echo "getting gff file for ${lcb}, chunk=${CHUNK} failed"; exit 1; }

                bedtools getfasta -s \
                    -bed "gff_chunks_${CHUNK}/${lcb}.gff" -fi all.fasta | \
                    sed 's/:.\+//g' > \
                    "fasta_chunks_${CHUNK}/${lcb}.fasta" \
                    || { echo "getting fasta file for ${lcb}, chunk=${CHUNK} failed"; exit 1; }
            done < "lcbs_to_align.${CHUNK}.txt"

            ## Align and produce chunk output file
            ## Another step required to put all output files together into an MAF
            "$maf_aligner" -l 1 -r 1 -e -8 \
                "gff_chunks_${CHUNK}" \
                "fasta_chunks_${CHUNK}" \
                "all.fasta.fai" \
                "lcbs_to_align.${CHUNK}.txt" \
                > "alignments/${CHUNK}.maf" \
                || { echo "aligning ${CHUNK} failed"; exit 1; }

            rm "lcbs_to_align.${CHUNK}.txt"
            rm -r "fasta_chunks_${CHUNK}" "gff_chunks_${CHUNK}"
            ;;

        "concatenate-maf")
#            cat alignments/*.maf > all.maf
#            gzip -f all.maf
#            rm all.fasta

            cp "$lcb_len" .
            "./$(basename "$lcb_len")" all.maf.gz > lcb_lengths.tsv \
                || { echo "getting LCB lengths failed"; exit 1; }
            ;;

        "variants-from-maf")
            cp "$varfinder" .

            ls "fasta" | sort | sed 's/.fasta//g' > "strains.txt" \
                || { echo "making list of strains failed"; exit 1; }

            "./$(basename "$varfinder")" "all.maf.gz" "strains.txt" > \
                "variants.tsv" \
                || { echo "finding variants failed"; exit 1; }
            gzip -f variants.tsv
            ;;

        "lcb-order-comparison")
            cp "$make_order_files" .
            cp "$compare_order" .

            ## Make commands
            Rscript -e 'x=scan("strains.txt",what="char");cc=combn(x,2);for(i in 1:ncol(cc)) cat(cc[1,i],"\t",cc[2,i],"\n",sep="");' \
                > "combinations.txt" \
                || { echo "making list of all strain combinations failed"; exit 1; }
            co="./$(basename "$compare_order")"
            rm -f commands
            while read -r strains; do
                s1="$(echo "$strains" | cut -f 1)"
                s2="$(echo "$strains" | cut -f 2)"
                oc="${co} \"order_files/${s1}.tsv\" \"order_files/${s2}.tsv\""
                echo "echo \"${s1}"$'\t'"${s2}"$'\t''$('"$oc"')"' >> \
                    commands \
                    || { echo "adding ${s1} and ${s2} to commands failed"; exit 1; }
            done < combinations.txt
            
            ## Make input files
            rm -rf order_files
            mkdir -p order_files
            zcat "lcb_table.size_filtered.tsv.gz" | \
                tail -n +2 | \
                cut -f 1 > \
                "filtered_lcbs.tsv" \
                || { echo "making list of size filtered lcbs failed"; exit 1; }
            gzip -f filtered_lcbs.tsv
            "./$(basename "$make_order_files")" \
                --output "order_files" \
                --lcbs "filtered_lcbs.tsv.gz" \
                "output.${KMER}.${MAX_BUBBLE_SIZE}/all_blocks.gff" \
                "strains.txt" \
                || { echo "making LCB order files failed"; exit 1; }

            ## Compare order for every pair
            cat "commands" | parallel -j$NTHREADS \
                > "order_comparison.txt" \
                || { echo "running order comparison in parallel failed"; exit 1; }

            ## Remove intermediate files
            rm -r combinations.txt commands order_files
            ;;

        "coord-table")
            ## UPDATE 2020-12-10: Added strand to output
            mkdir -p coords
            "$coord_table" all.maf.gz "$CHUNK" \
                | sort -t$'\t' -k 1n,1 -k 2n,2 -k 8rn,8 \
                > "coords/${CHUNK}.tsv" \
                || { echo "converting alignments for ${CHUNK} to table failed"; exit 1; }
            gzip --best -f "coords/${CHUNK}.tsv"
            ;;

        "overlap-orthologs-1")
            ## Use coords table results to infer orthologs
            ## Part 1: per strain, find coordinates for each strain
            mkdir -p overlap
            overlapdir="${SCRATCHDIR}/${SPECIES}-overlap"
            mkdir -p "$overlapdir"

            ## Input files
            coords="coords/${CHUNK/__pbnb/__pbmg}.tsv.gz"
            gff="gff/${CHUNK}.gff3"

            ## Modify gff file for the NCBI PacBio annotations
            if [[ "$gff" =~ __pbnb ]]; then
                ## Change chromosome names
                echo "CWD $(pwd)"
                python -c "$gff_filter" "$gff"  \
                    | python -c "$change_replicon_names" "gff/${CHUNK/pbnb/pbmg}.gff3" \
                    | sed 's/.\+__.\{2,7\}\.//g' \
                    | awk -F$'\t' '$3 == "gene" { print };' \
                    > "${overlapdir}/${CHUNK}.gff3" \
                    || { echo "filtering to genes for ${CHUNK} failed"; exit 1; }
            else
                python -c "$gff_filter" "$gff"  \
                    | sed 's/.\+__.\{2,7\}\.//g' \
                    | sed 's/^NODE/lcl|NODE/g' \
                    | awk -F$'\t' '$3 == "gene" { print };' \
                    > "${overlapdir}/${CHUNK}.gff3" \
                    || { echo "filtering to genes for ${CHUNK} failed"; exit 1; }
            fi

            ## Find gene locations
            zcat "$coords" \
                | awk -F$'\t' '$8 == 1 {s="+"}; $8 == -1 {s="-"}; $6 != $7 {print $5 "\t" $6 "\t" $7 "\t" $1 "_" $2 "_" $3 "\t0\t" s}' \
                | bedtools sort -i - \
                > "${overlapdir}/${CHUNK}.input.bed" \
                || { echo "getting LCB locations for ${CHUNK} failed"; exit 1; }

            bedtools sort -i "${overlapdir}/${CHUNK}.gff3" \
                > "${overlapdir}/${CHUNK}.sorted.gff3" \
                || { echo "sorting gff for ${CHUNK} failed"; exit 1; }

            ## TODO: "inconsistent naming convention" is missing genes
            ## 2021-06-04: Is this still true?
            bedtools intersect -wb -sorted -nonamecheck \
                -a "${overlapdir}/${CHUNK}.input.bed" \
                -b "${overlapdir}/${CHUNK}.sorted.gff3" \
                | cut -f 4,6,10,11,13,15 \
                | sort -k 1b,1 -T "$SCRATCHDIR" \
                | awk '{print $1 "\t" $2 "\t" $5 "\t" $6 "\t" "'"$CHUNK"'" "\t" $4-$3+1};' \
                | python -c "$extract_locus_tag" \
                > "${overlapdir}/annotated.${CHUNK}.tsv" \
                || { echo "finding gene locations for ${CHUNK} failed"; exit 1; }

            rm "${overlapdir}/$CHUNK".*
            ;;

        "overlap-orthologs-2")
            ## Part 2: per strain count overlaps with other strains
            overlapdir="${SCRATCHDIR}/${SPECIES}-overlap"
            "$tally_overlap" \
                "${overlapdir}/annotated.${CHUNK}.tsv" \
                ${overlapdir}/annotated.*.tsv \
                > "${overlapdir}/${CHUNK}.overlap_tally.tsv" \
                || { echo "Counting overlaps for ${CHUNK} failed"; exit 1; }
            ;;

        "overlap-orthologs-3")
            ## First, filter to just the connections that meet the
            ## threshold
            overlapdir="${SCRATCHDIR}/${SPECIES}-overlap"
            THRESHOLD=0.80
#            cp "$overlap_tally_to_ologs" .
#            "./$(basename "$overlap_tally_to_ologs")" $THRESHOLD \
#                "$overlapdir"/*.overlap_tally.tsv \
#                > "${overlapdir}/graph.ncol" \
#                || { echo "getting graph input for ${SPECIES} failed"; exit 1; }
#
#            ## Then, print out connnected components; uses igraph.
#            ## Also adds in genes that are not in the graph
#            ls "$overlapdir" | grep ".overlap_tally.tsv" \
#                | sed 's/\.overlap_tally.tsv//g' \
#                > assemblies.txt
#            ## Updated 2021-06-03 with code that does not add in a bunch
#            ## of fragmented singleton genes at the end. The output file
#            ## of the updated was named ".filtered" so that it would not
#            ## interfere with existing scripts.
#            ## TODO: If this script is cleaned up and re-run, the
#            ## output files of the next two steps do not need to have
#            ## .filtered in their names and the commented out parts can
#            ## be removed.
#            #python -c "$print_orthologs" "${overlapdir}/graph.ncol" \
#            #    assemblies.txt gff \
#            python -c "$print_orthologs" "${overlapdir}/graph.ncol" \
#                assemblies.txt \
#                > "overlap_orthologs.${THRESHOLD}.filtered.tsv" \
#                || { echo "making ortholog table for ${SPECIES} failed"; exit 1; }
#                #> "overlap_orthologs.${THRESHOLD}.tsv" \
#                #|| { echo "making ortholog table for ${SPECIES} failed"; exit 1; }

            ## Updated 2021-04-01 to identify more fragmented genes and pseudogenes
            ## However, the pseudogenes and fragments should have been filtered
            ## out at the beginning, so this is just an extra check.
            ## Updated again 2021-06-03: this was marking some genes as partial
            ## that are not. Fixed (probably). Output was named .filtered
            ## to avoid interfering with existing pipelines.
            ## Updated 2021-06-15 to include start and end positions.
            ## Re-run 2021-07-30 on actually filtered data!
            echo "$annotation_table_from_ologs" > annot_table_ologs.py
            chmod u+x annot_table_ologs.py
            ./annot_table_ologs.py \
                "overlap_orthologs.${THRESHOLD}.filtered.tsv" \
                annotation.tsv.gz \
                | sort -t$'\t' -k 1b,1 -k 2b,2 \
                > "annotation_table.overlap.${THRESHOLD}.filtered.tsv" \
                || { echo "annotating overlap orthologs failed"; exit 1; }
                #> "annotation_table.overlap.${THRESHOLD}.tsv" \
                #|| { echo "annotating overlap orthologs failed"; exit 1; }
            gzip -f "annotation_table.overlap.${THRESHOLD}.filtered.tsv" \
                || { echo "gzipping overlap orthologs annotation failed"; exit 1; }
            #gzip -f "annotation_table.overlap.${THRESHOLD}.tsv" \
            #    || { echo "gzipping overlap orthologs annotation failed"; exit 1; }
#            cp "${overlapdir}/graph.ncol" . # Is it too big? Move to second tier
            ;;

        "make-tree")
            cd "$OUTDIR"
            rm -rf _
            mkdir -p parsnp-tree && cd parsnp-tree
            echo "PWD: $(pwd)"

            ## Copy sequences of genomes. Keep only one assembly per strain,
            ## with priority on our version of the assembly and on more
            ## complete assemblies.
            rm -rf fasta
            mkdir -p fasta
            for species in "${SPECIES_LIST[@]}"; do
                for suffix in __pbmg __inr __inb __magncbi __crfs; do
                    for fname in "../${species}/fasta/"*"$suffix".fasta; do
                        echo "$fname"
                        strain="$(basename "$fname" | sed 's/'"$suffix"'.fasta//g')"
                        if [[ -f "fasta/${strain}.fasta" ]]; then continue; fi
                        cp "$fname" "fasta/${strain}.fasta" \
                            || { echo "copying ${strain} from ${fname} failed"; exit 1; }
                    done
                done
            done

            ## Run with v1.2
#            parsnp -r fasta/1021.fasta -d "./fasta" \
#                -p "$NTHREADS" -x -D 0.16 -c -o "output" \
#                || { echo "running parsnp to make tree failed"; exit 1; }

            mv fasta/1021.fasta reference.fasta
            for replicon in chromosome pSymA pSymB; do
                blgrep "$replicon" reference.fasta \
                    > fasta/1021.fasta \
                    || { echo "extracting ${replicon} failed"; exit 1; }
                ## Run with v1.5.6
                parsnp --use-fasttree -r fasta/1021.fasta -d "./fasta" \
                    -p "$NTHREADS" -x -D 0.16 -c -o "output-${replicon}" \
                    || { echo "running parsnp to make tree for ${replicon} failed"; exit 1; }
            done

            rm -rf fasta reference.fasta
            ;;

    esac

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes ($(($SECONDS/3600))) hours)"
    echo "done"

    rm "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-${1}-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    case "$1" in

        "get-sequences")
            WALLTIME="00:10:00"
            MEM="1GB"
            NTHREADS=1
            for species in "${SPECIES_LIST[@]}"; do
                for subdir in "${!STRAIN_PREFIXES[@]}"; do
                    echo "$1"$'\t'"$species"$'\t'"$subdir" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "align")
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'"_" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "process")
            WALLTIME=02:00:00
            NTHREADS=3
            MEM=16GB
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'"_" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "maf")
            WALLTIME=04:00:00
            NTHREADS=1
            MEM=4GB
            for species in "${SPECIES_LIST[@]}"; do
                for chunk in $(seq 1 $N_CHUNKS); do
                    echo "$1"$'\t'"$species"$'\t'"_"$'\t'"$chunk" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "concatenate-maf")
            WALLTIME=00:05:00
            NTHREADS=1
            MEM=4GB
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'"_"$'\t'"_" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "variants-from-maf")
            WALLTIME=01:00:00
            NTHREADS=1
            MEM=4GB
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'"_"$'\t'"_" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "lcb-order-comparison")
            WALLTIME=03:00:00
            MEM=62GB
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'"_"$'\t'"_" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "coord-table")
            WALLTIME=00:10:00
            NTHREADS=1
            MEM=1GB
            for species in "${SPECIES_LIST[@]}"; do
                for ffile in "${OUTDIR}/${species}/fasta/"*; do
                    strain="$(basename "${ffile/.fasta/}")"
                    echo "$1"$'\t'"$species"$'\t'"_"$'\t'"$strain" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "overlap-orthologs-1")
            WALLTIME=02:00:00 ## Nearly all will finish in < 10 minutes
            NTHREADS=1
            MEM=4GB
            for species in "${SPECIES_LIST[@]}"; do
                for gfile in "${OUTDIR}/${species}/gff/"*; do
                    strain="$(basename "${gfile/.gff3/}")"
                    if [[ "$strain" =~ rast ]]; then continue; fi
                    echo "$1"$'\t'"$species"$'\t'"_"$'\t'"$strain" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "overlap-orthologs-2")
            WALLTIME=08:00:00
            NTHREADS=1
            MEM=16GB
            for species in "${SPECIES_LIST[@]}"; do
                for gfile in "${OUTDIR}/${species}/gff/"*; do
                    strain="$(basename "${gfile/.gff3/}")"
                    if [[ "$strain" =~ rast ]]; then continue; fi
                    echo "$1"$'\t'"$species"$'\t'"_"$'\t'"$strain" > "${ARRAYDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

        "overlap-orthologs-3")
            WALLTIME=03:00:00
            NTHREADS=1
            MEM=64GB
            for species in "${SPECIES_LIST[@]}"; do
                echo "$1"$'\t'"$species"$'\t'_$'\t'_ > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        "make-tree")
            WALLTIME="06:00:00"
            MEM=64GB
            echo "$1"$'\t'_$'\t'_$'\t'_ > "${ARRAYDIR}/${i}"
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
        --partition="$QUEUE" --time=$WALLTIME --mem=$MEM \
        --nodes=1 --ntasks=$NTHREADS --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
