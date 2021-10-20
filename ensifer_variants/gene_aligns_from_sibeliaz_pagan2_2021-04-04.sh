#!/bin/bash
#
# Create alignments of orthologous genes identified by an overlap cutoff
# from SibeliaZ alignments.
#
# Protein-coding genes are aligned by using a codon alignment, while
# other genes are aligned by their nucleotide sequences.
#
# Essentially identical to the 2021-02-10 run, but I made a new script
# because I updated the orthologs, and I wanted to stay organized.
#
# TODO: The input for this script was updated 2021-06-03 to better handle
# fragment filtering:
#   - singleton genes that are fragments were removed.
#   - a few genes labeled as fragments that aren't (NCBI PacBio strains)
#     are no longer labeled as fragments.
#   - At the time, this had little to no effect on any downstream analyses,
#     but this script should eventually be replaced with one that takes
#     in the corrected input data.
# NOTE: The 2020-06-07 uses the corrected input. Differences are quite small.
#
# SETTINGS
#
QUEUE=amdsmall
#
N_BATCHES=150
N_BATCHES_NUC=8
SEED=681
NTHREADS_BLAST=32
MIN_ID=80

PROJDIR="${HOME}/project/jgi_sequencing"
PBNB_ASSEMBLY_DIR="${PROJDIR}/data/assembly/complete/ncbi_pacbio_2016"
ORTHODIR="${PROJDIR}/results/ortho/medmel/sibeliaz/2020-04-09_split"
ORTHOLOG_FILE=overlap_orthologs.0.80.tsv
ORTHOLOG_ANNOTATION=annotation_table.overlap.0.80.tsv.gz
GENOME_DIR=fasta
GFF_DIR=gff

OUTDIR="${PROJDIR}/results/gene_aligns/from_sibeliaz/pagan2/2021-04-04"
ARRAYDIR="${OUTDIR}/arrayjobdata"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

BINDIR="${PROJDIR}/script/bin"
cons="${BINDIR}/fasta_consensus"
fasta2vcf="${BINDIR}/fasta_to_vcf"
extract_hits="${BINDIR}/extract_blast_fmt7_lengths_bitscore.py"
bbh="${BINDIR}/extract_bbh.py"

filter_singletons="
## Filter out singleton sequences
## Could be improved by keeping track of different annotations
## and assemblies of the same strain.
import csv
import sys
import pyfastx
keep = set()
with open(sys.argv[3], 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        n = 0
        for c, v in row.items():
            if c == 'ortholog': continue
            if v != '': n += len(v.split(','))
        if n > 1: keep.add(row['ortholog'])
for name, seq, _ in pyfastx.Fastx(sys.argv[2]):
    if name in keep:
        print('>' + sys.argv[1] + '::' + name)
        print(seq)
"

change_replicon_names="
import sys
from Bio import SeqIO
old_sizes = {}
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    old_sizes[len(rec)] = rec.id
for rec in SeqIO.parse(sys.argv[2], 'fasta'):
    l = len(rec)
    for offset in range(-10, 10, 1):
        if (l + offset) in old_sizes:
            print('>' + old_sizes[l+offset])
            print(str(rec.seq))
            break
"

prep_gff="#!/usr/bin/env python3
import sys
locus_tags = {}
biotypes = {}
keep = set()
with open(sys.argv[3], 'rt') as handle:
    for line in handle:
        if line[0] == '#': continue
        fields = line.strip().split('\t')
        tags = {x.split('=',1)[0]:x.split('=',1)[1] for x in fields[8].split(';')}
        key = (fields[0], fields[3], fields[4], fields[6])
        if fields[2] == 'gene': keep.add(key)
        if 'locus_tag' in tags: locus_tags[key] = tags['locus_tag']
        if 'gene_biotype' in tags:
            if tags['gene_biotype'] == 'protein_coding': biotypes[key] = 'p'
            else: biotypes[key] = 'n'
        if fields[2] == 'CDS': biotypes[key] = 'p'
for key in keep:
    name = sys.argv[1] + '_' + locus_tags[key]
    if sys.argv[2] == 'n' and (key in biotypes and biotypes[key] != 'n'): continue
    if sys.argv[2] == 'p' and (key not in biotypes or biotypes[key] == 'n'): continue
    print(key[0] + '\t' + str(int(key[1])-1) + '\t' + key[2] + '\t' + name + '\t.\t' + key[3])
"

## bedtools fasta doesn't handle circular sequences, so I had to make
## my own version.
get_fasta="#!/usr/bin/env python3
import sys
from Bio import SeqIO, Seq
seqs = {}
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    seqs[rec.id] = str(rec.seq)
for line in sys.stdin:
    chrom, s, e, name, _, strand = line.strip().split('\t')
    l = len(seqs[chrom])
    if int(e) > l:
        s = seqs[chrom][int(s):] + seqs[chrom][0:(int(e)-l)]
    else:
        s = seqs[chrom][int(s):int(e)]
    if strand == '-':
        s = Seq.reverse_complement(s)
    print('>' + name + '\n' + s)
"

make_unaligned_files="#!/usr/bin/env python3
import csv
import sys
from Bio import SeqIO
prot = SeqIO.index('protein.fasta', 'fasta')
nuc = SeqIO.index('nucleotide.fasta', 'fasta')
exclude = set()
with open('partials_and_pseudogenes.txt', 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        exclude.add(tuple(row))
with open(sys.argv[1], 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        ol = row['ortholog']
        sys.stderr.write('RUNNING %s\n' % ol)
        p = False
        n = False
        nucs = set(); prots = set()
        s = []
        for strain, genes in row.items():
            st = strain.split('__')[0]
            if strain == 'ortholog': continue
            for gene in genes.split(','):
                name = strain + '_' + gene
                if (ol, st, gene) in exclude:
                    continue
                if name in prot:
                    p = True
                    s.append((name, str(prot[name].seq)))
                    prots.add(name)
                elif name in nuc:
                    n = True
                    s.append((name, str(nuc[name].seq)))
                    nucs.add(name)
        if p and n:
            ## This seems to happen when the gene is split into multiple
            ## pieces such that the CDS and gene entry coordinates are
            ## not identical. Mostly transposases with ribosomal slippage.
            sys.stderr.write('SKIPPING: %s, unclear whether protein or nucleotide\n' % ol)
            sys.stderr.write('\tNUCLEOTIDE SEQS: ' + '\t'.join(nucs) + '\n')
            sys.stderr.write('\tPROTEIN SEQS: ' + '\t'.join(prots) + '\n')
        else:
            if p: fname = 'unaligned-protein/' + ol + '.fasta'
            elif n: fname = 'unaligned-nucleotide/' + ol + '.fasta'
            if len(s) == 0:
                sys.stderr.write('Empty list of sequences for %s\n' % ol)
            sys.stderr.write('PRINTING %s into %s\n' % (ol, fname))
            with open(fname, 'wt') as ohandle:
                for name, sq in s:
                    if p: sq += 'N' * (len(sq) % 3)
                    ohandle.write('>' + name + '\n' + sq.upper() + '\n')
"

trim_terminal_n="
import sys
from Bio import AlignIO
a = AlignIO.read(sys.argv[1], 'fasta')
if len(a) == 1: exit()
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
with open(sys.argv[1], 'wt') as out:
    for rec in a:
        if rec.id.endswith('/1'):
            rec.id = rec.id[:-2]
        elif rec.id.endswith('/2') or rec.id.endswith('/3'):
            continue
        out.write('>' + rec.id + '\n' + str(rec.seq) + '\n')
"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load pagan2/20190829
    module load mafft/7.475
    module load fasttree/2.1.8
    module load bowtie2/2.4.2
    module load samtools/1.10
    module load mash/2.1
    module load gcc
    module load bcftools/1.10.2
    module load ncbi_blast+/2.8.1

    set -euo pipefail

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    SPECIES="$(cut -f 2 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    BATCH="$(cut -f 3 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    cd "$OUTDIR"

    mkdir -p "$SPECIES" && cd "$SPECIES"

    case "$TASK" in

        "prep")
            ## Clear directories to hold sequences
#            rm -rf unaligned-protein unaligned-nucleotide aligned-protein aligned-nucleotide
            mkdir -p unaligned-protein unaligned-nucleotide aligned-protein aligned-nucleotide

            orthologs="${ORTHODIR}/${SPECIES}/${ORTHOLOG_FILE}"

#            ## Make a master file of protein-coding and non-protein-coding 
#            ## gene sequences
#            echo "$prep_gff" > prep_gff.py; chmod u+x prep_gff.py
#            echo "$get_fasta" > get_fasta.py; chmod u+x get_fasta.py
#            rm -f nucleotide.fasta protein.fasta
#            for strain in $(awk 'NR == 1 { for(i=2;i<=NF;i++) { print $i }; };' "$orthologs"); do
#                gff="${ORTHODIR}/${SPECIES}/${GFF_DIR}/${strain}.gff3"
#                if [[ "$strain" =~ __pbnb ]]; then
#                    ## MaGe and NCBI assemblies can be very slightly different; Get NCBI
#                    ## assembly here and rename replicons to match the mage assembly
#                    mage_fasta="${ORTHODIR}/${SPECIES}/${GENOME_DIR}/${strain/__pbnb/__pbmg}.fasta"
#                    ncbi_fasta="$(echo "${PBNB_ASSEMBLY_DIR}/${strain/_1/-1}/genome.fasta" | sed 's/__pbnb//g')"
#                    python3 -c "$change_replicon_names" "$mage_fasta" "$ncbi_fasta" \
#                        > tmp.fasta \
#                        || { echo "changing replicon names for ${strain} failed"; exit 1; }
#                    fasta=tmp.fasta
#                else
#                    fasta="${ORTHODIR}/${SPECIES}/${GENOME_DIR}/${strain}.fasta"
#                fi
#                ./prep_gff.py "$strain" n "$gff" \
#                    | ./get_fasta.py "$fasta" \
#                    >> nucleotide.fasta \
#                    || { echo "adding nucleotide gene sequences from ${strain} failed"; exit 1; }
#                ./prep_gff.py "$strain" p "$gff" \
#                    | ./get_fasta.py "$fasta" \
#                    >> protein.fasta \
#                    || { echo "adding protein gene sequences from ${strain} failed"; exit 1; }
#            done
#
#            annot="${ORTHODIR}/${SPECIES}/${ORTHOLOG_ANNOTATION}"
            ## TODO if re-run: there are fragments/partials among the singleton
            ## genes (though not all of them), but this table also mistakenly
            ## labels a few genes as partial that are not.
#            zcat "$annot" \
#                | awk -F$'\t' '$8 == 1 || $7 == 1 { print $1 "\t" $2 "\t" $5 };' \
#                > partials_and_pseudogenes.txt \
#                || { echo "making list of genes to exclude failed"; exit 1; }

            ## For each ortholog, make a file of unaligned sequences
            echo "$make_unaligned_files" > make_files.py; chmod u+x make_files.py
            ./make_files.py "$orthologs" \
                || { echo "making unaligned files for ${SPECIES} failed"; exit 1; }

            ;;

        "align-protein-coding")
            ls unaligned-protein | awk '{ print "unaligned-protein/" $1};' > "files.p.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES" "files.p.${BATCH}.txt" > \
                "files.only.p.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} in ${GROUP} failed"; exit 1; }

            ## A few notes:
            ## 1. Have to load the mafft and fasttree modules after pagan2
            ##    so the right version of those programs gets used
            ## 2. For reasons that aren't clear, pagan2 cannot run its
            ##    own mafft or fasttree versions during a batch script
            ##    (except very occasionally, it can), so these have to
            ##    be run separately
            ## 3. The issue has something to do with /proc/self/exe
            ##    not pointing to the right spot, probably. Too annoying
            ##    to figure out.
            while read -r fname; do
                echo "RUNNING ${fname}"
                outprefix="aligned-protein/$(basename "${fname/.fasta/}")"
                if [[ -f "${outprefix}.fas" && -s "${outprefix}.fas" ]]; then
                    echo "already finished ${fname} in ${SPECIES}, ${BATCH}"
                    continue
                fi
                if [[ "$(grep -c "^>" "$fname")" == "1" ]]; then
                    cp "$fname" "${outprefix}.fas"
                    echo "only one sequence in ${fname}, copying and skipping"
                    continue
                fi
                grep "^>" "$fname" > /dev/null || { echo "${fname} empty, skipping"; continue; }
                mafft --auto --quiet "$fname" > "tmp.${BATCH}" \
                    || { echo "aligning ${fname} with MAFFT failed"; exit 1; }
                fasttree -quiet -nt "tmp.${BATCH}" > "${outprefix}.tre" \
                    || { echo "making tree for ${fname} with FastTree failed"; exit 1; }
                Rscript -e 'library(ape); write.tree(multi2di(read.tree("'"${outprefix}.tre"'")), "'"tmp.${BATCH}.tre"'");' \
                    || { echo "converting polytomies to dichotomies for ${fname} in ${SPECIES}, ${BATCH} failed"; exit 1; }
                pagan2 --silent --codons -o "$outprefix" -s "$fname" -t "tmp.${BATCH}.tre" \
                    || { echo "alignment on ${fname} for ${SPECIES} (batch ${BATCH}) failed"; exit 1; }
                python -c "$trim_terminal_n" "${outprefix}.fas" \
                    || { echo "trimming ${fname} alignment for ${SPECIES} (batch ${BATCH}) failed"; exit 1; }
            done < "files.only.p.${BATCH}.txt"

            ## Delete the unaligned sequences (do this manually)
            #rm -r "unaligned-aa"
            rm -f "files.p.${BATCH}.txt" "files.only.p.${BATCH}.txt" "tmp.${BATCH}" "tmp.${BATCH}.tre"
            ;;

        "align-nucleotide")
            ls unaligned-nucleotide | awk '{ print "unaligned-nucleotide/" $1};' > "files.n.${BATCH}.txt"
            split -n "r/$BATCH/$N_BATCHES_NUC" "files.n.${BATCH}.txt" > \
                "files.only.n.${BATCH}.txt" \
                || { echo "splitting files for ${BATCH} in ${GROUP} failed"; exit 1; }

            while read -r fname; do
                echo "RUNNING ${fname}"
                outprefix="aligned-nucleotide/$(basename "${fname/.fasta/}")"
                wrong_outprefix="aligned-protein/$(basename "${fname/.fasta/}")"
                rm -f "${wrong_outprefix}.fas" "${wrong_outprefix}.tre"
                if [[ -f "${outprefix}.fas" && -s "${outprefix}.fas" ]]; then
                    echo "already finished ${fname} in ${SPECIES}, ${BATCH}"
                    continue
                fi
                if [[ "$(grep -c "^>" "$fname")" == "1" ]]; then
                    cp "$fname" "${outprefix}.fas"
                    echo "only one sequence in ${fname}, copying and skipping"
                    continue
                fi
                grep "^>" "$fname" > /dev/null || { echo "${fname} empty, skipping"; continue; }
                mafft --auto --quiet "$fname" > "tmp.${BATCH}" \
                    || { echo "aligning ${fname} with MAFFT failed"; exit 1; }
                fasttree -quiet -nt "tmp.${BATCH}" > "${outprefix}.tre" \
                    || { echo "making tree for ${fname} with FastTree failed"; exit 1; }
                Rscript -e 'library(ape); write.tree(multi2di(read.tree("'"${outprefix}.tre"'")), "'"tmp.${BATCH}.tre"'");' \
                    || { echo "converting polytomies to dichotomies for ${fname} in ${SPECIES}, ${BATCH} failed"; exit 1; }
                pagan2 --silent -o "$outprefix" -s "$fname" -t "tmp.${BATCH}.tre" \
                    || { echo "alignment on ${fname} for ${SPECIES} (batch ${BATCH}) failed"; exit 1; }
                python -c "$trim_terminal_n" "${outprefix}.fas" \
                    || { echo "trimming ${fname} alignment for ${SPECIES} (batch ${BATCH}) failed"; exit 1; }
            done < "files.only.n.${BATCH}.txt"

            ## Delete the unaligned sequences (do this manually)
            #rm -r "unaligned-nuc"
            rm -f "files.n.${BATCH}.txt" "files.only.n.${BATCH}.txt" "tmp.${BATCH}" "tmp.${BATCH}.tre"

            ;;

        consensus)
            rm -f consensus.fasta
            for fname in aligned-nucleotide/*.fas aligned-protein/*.fas; do
                if [[ ! -s "$fname" ]]; then continue; fi
                "$cons" "$fname" "$(basename "${fname/.fas/}")" $SEED \
                    >> consensus.fasta \
                    || { echo "getting consensus sequence for ${fname} failed"; exit 1; }
            done

            rm -rf vcf
            mkdir -p vcf
            ## NOTE: This reports an allele for every gene sequence;
            ## not organized by strain.
            for fname in aligned-protein/*.fas aligned-nucleotide/*.fas; do
                if [[ ! -s "$fname" ]]; then continue; fi
                "$fasta2vcf" "$fname" > "vcf/$(basename "${fname/.fas/.vcf}")" \
                    || { echo "running fasta2vcf on ${fname} in ${SPECIES} failed"; exit 1; }
            done
            ;;

        vcf-spanfran)
            ## Filter to just the SpanFran strains and get rid of non-segregating
            ## sites
            rm -rf spanfran_vcf
            mkdir -p spanfran_vcf
            for fname in aligned-protein/*.fas aligned-nucleotide/*.fas; do
                if [[ ! -s "$fname" ]]; then continue; fi
                blgrep "__magncbi" "$fname" > tmp || continue
                if [[ ! -s tmp ]]; then continue; fi
                "$fasta2vcf" tmp \
                    | bcftools view -Ov -m2 --min-ac 1:nonmajor \
                    > "spanfran_vcf/$(basename "${fname/.fas/.vcf}")" \
                    || { echo "running fasta2vcf on ${fname} in ${SPECIES} failed"; exit 1; }
            done
            rm tmp
            ;;

        species-orthologs-bbh)
            ## Check for similarity among all non-singleton genes
            ## using consensus sequences

            ## Filter out genes with just one sequence, because I think
            ## these are usually fragments
            for species in meliloti medicae; do
                python -c "$filter_singletons" \
                    "$species" \
                    "../${species}/consensus.fasta" \
                    "${ORTHODIR}/${species}/${ORTHOLOG_FILE}" \
                    > "${species}.fasta" \
                    || { echo "assembling sequences from ${species} failed"; exit 1; }
                ## Make blast database
                makeblastdb -dbtype nucl -in "${species}.fasta" -out "${species}.cons" \
                    || { echo "making blast DB failed"; exit 1; }
            done

            ## Run blast in each direction
            rm -f blast_results.txt
            blastn -query meliloti.fasta -db medicae.cons \
                -outfmt '7 std qlen slen' \
                -num_threads $NTHREADS_BLAST \
                -soft_masking true \
                -evalue 1e-7 \
                >> blast_results.txt \
                || { echo "blasting meliloti against medicae failed"; exit 1; }
            blastn -query medicae.fasta -db meliloti.cons \
                -outfmt '7 std qlen slen' \
                -num_threads $NTHREADS_BLAST \
                -soft_masking true \
                -evalue 1e-7 \
                >> blast_results.txt \
                || { echo "blasting medicae against meliloti failed"; exit 1; }

            ## Apply a minimum % identity filter to prevent nonsense
            rm -f all_v_all.filtered.txt
            awk 'BEGIN { FS="\t" }; $3 >= '$MIN_ID' {print};' blast_results.txt \
                > filtered.txt \
                || { echo "applying minimum %id filter failed"; exit 1; }

            cp "$extract_hits" .
            cp "$bbh" .
            # Extract hits
            ./$(basename "$extract_hits") filtered.txt \
                > hits.tsv \
                || { echo "extracting hits failed"; exit 1; }
            # Get BBHs
            ./$(basename "$bbh") hits.tsv \
                > bbh.tsv \
                || { echo "getting BBHs failed"; exit 1; }

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
            WALLTIME=01:30:00
            MEM=2GB
            NTHREADS=1
            for species in medicae meliloti; do
                echo "$1"$'\t'"$species"$'\t'_ > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        align-*)
            WALLTIME=05:00:00
            MEM=8GB
            NTHREADS=1
            case "$1" in
                align-nucleotide) N_BATCHES=$N_BATCHES_NUC ;;
            *) N_BATCHES=$N_BATCHES ;;
        esac
        for species in medicae meliloti; do
            for batch in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'"$species"$'\t'"$batch" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
        done
        ;;

    consensus)
        WALLTIME=06:00:00
        MEM=8GB
        NTHREADS=1
        for species in medicae meliloti; do
            echo "$1"$'\t'"$species"$'\t'_ > "${ARRAYDIR}/${i}"
            i=$(($i+1))
        done
        ;;

    vcf-spanfran)
        WALLTIME=01:00:00
        MEM=8GB
        NTHREADS=1
        for species in medicae meliloti; do
            echo "$1"$'\t'"$species"$'\t'_ > "${ARRAYDIR}/${i}"
            i=$(($i+1))
        done
        ;;

    species-orthologs-bbh)
        WALLTIME=00:25:00
        MEM=16GB
        NTHREADS=$NTHREADS_BLAST
        echo "$1"$'\t'bbh$'\t'_ > "${ARRAYDIR}/${i}"
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
    --partition=$QUEUE --nodes=1 --ntasks=$NTHREADS \
    --time=$WALLTIME --mem=$MEM --array="$array" "$sfile"
echo "$WORKDIR"
echo "Submitted ${i} jobs"

fi
