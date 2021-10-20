#!/bin/bash
#
# Run allele-specific fitness on full dataset.
#
# RDVs are filtered out.
#
# UPDATED 2021-08-26: New v1.8 annotation used for gene annotations (not for
# variant filtering)
#
# SETTINGS
#
# pbs
QUEUE="amdsmall"
#
MIN_HOST_MAF=0.1
THINNING_WINDOW="2 kb"
THINNING_STEP=20
THINNING_THRESHOLD=0.90
MIN_SYM_MAF=0.1 # Applied to nodule communities
MAX_MISS=0.2
ENS_R2_THRESHOLD=0.95
N_BATCHES=500
GENE_WINDOW_FILTER=0 # Window around Mtr genes in which to include SNPs in the analysis
MIN_FOLD_CHANGE=1.5  # Minimum fold change in fitness to count as a candidate
RANDOM_FRACTION=4E-06 # Fraction of pairs to randomly keep
LD_GROUP_SIZE=3

combine_partners="#!/usr/bin/env python3
import csv
import sys

mtr = {}
with open(sys.argv[1], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        mtr[row[2]] = (row[0], row[5], row[7])

print('ensifer_variant', 'medicago_variant', 'ld_group',
      'rank', 'log2_fitness_diff', 'ensifer_gene',
      'usda1106_ncbi', 'usda1106_mage', 'ensifer_annotations',
      'medicago_gene', 'medicago_annotations', sep='\t')
with open(sys.argv[2], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        rank = row[2]
        print(row[0], mtr[rank][0], row[1], rank,
              row[3], row[4], row[5], row[6], row[7],
              mtr[rank][1], mtr[rank][2], sep='\t')
"

combine_ensifer_annotations="
import csv
import itertools
import sys
for rs, rows in itertools.groupby(csv.reader(sys.stdin, delimiter='\t'), lambda r: (r[0], r[2])):
    prods = set()
    mage = set()
    ncbi = set()
    ols = set()
    for row in rows:
        if row[4] != '': ols.update(row[4].split(','))
        if row[5] != '': mage.update(row[5].split(','))
        if row[6] != '': ncbi.update(row[6].split(','))
        if row[7] != '': prods.update(row[7].split(';'))
    row[4] = ','.join(ols)
    row[5] = ','.join(mage)
    row[6] = ','.join(ncbi)
    row[7] = ';'.join(prods)
    print('\t'.join(row))
"

list_ld_groups="
import sys
groups = {}
with open(sys.argv[2], 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        groups[row[3]] = 'group-' + row[4]
with open(sys.argv[1], 'rt') as handle:
    for line in handle:
        print(groups[line.strip()])

"

add_ld_groups="
import sys
groups = {}
with open(sys.argv[2], 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        groups[row[3]] = ('group-' + row[4], set(row[5].split(',')))
with open(sys.argv[1], 'rt') as handle:
    for i, line in enumerate(handle):
        row = line.strip().split('\t')
        g = groups[row[0]]
        for rs in g[1]:
            print(rs, g[0], i+1, row[14], sep='\t')
"

annotate="#!/usr/bin/env python3
import csv
import gzip
import itertools
import sys
ologs = set()
annotation_info = []
with open(sys.argv[1],'rt') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        variant = row[0]
        olog = variant.split('-')[1]
        ologs.add(olog)
        annotation_info.append([olog] + row)
with gzip.open(sys.argv[2], 'rt') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for olog, rows in itertools.groupby(rdr, lambda r: r[0]):
        if olog not in ologs: continue
        usda1106_lts_ncbi = set()
        usda1106_lts_mage = set()
        products = set()
        for row in rows:
            if row[6] == '1': continue
            products.add(row[5])
            if row[1] == 'USDA1106':
                if 'CDO' in row[4]:
                    usda1106_lts_ncbi.add(row[4])
                elif 'SMEL' in row[4]:
                    usda1106_lts_mage.add(row[4])
        for x in annotation_info:
            if x[0] == olog:
                x += [','.join(usda1106_lts_ncbi), ','.join(usda1106_lts_mage), ';'.join(products)]
## variant ID, group, rank, fitdiff, OL, NCBI lt, MaGe lt, descriptions
for x in annotation_info:
    print(x[1], x[2], x[3], x[4], x[0], x[5], x[6], x[7], sep='\t')
"

format_medicago="
import collections
import sys
products = collections.defaultdict(set)
lts = collections.defaultdict(set)
info = {}
for line in sys.stdin:
    row = line.strip().split('\t')
    rank, rs, fd = row[3].split(':')
    if len(row) > 12:
        lt = row[12]
        if lt != '' and lt != '.': lts[rs].add(lt)
    info[int(rank)] = (rs, fd)
    if len(row) > 13:
        prod = row[13]
        if prod != '' and prod != '.': products[rs].add(prod)
for rank in sorted(info):
    rs = info[rank][0]
    print(rs, '', rank, info[rank][1], '',
          ','.join(lts[rs]), '', ';'.join(products[rs]), sep='\t')
"

filter_candidates="#!/usr/bin/env python3
import csv
import math
import random
import sys
threshold = math.log(${MIN_FOLD_CHANGE}, 2)
chance = float(${RANDOM_FRACTION})
rand = random.Random()
with open(sys.argv[1] + '.top.tsv', 'wt') as ch:
    with open(sys.argv[1] + '.random.tsv', 'wt') as rh:
        for fname in sys.argv[2:]:
            with open(fname, 'rt') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    if row['host_variant'].startswith('rdv'): continue
                    l2 = math.log(float(row['wAB']), 2) - math.log(float(row['wAb']), 2)
                    if abs(l2) >= threshold:
                        ch.write('\t'.join(row.values()) + '\t' + str(l2) +
                                 '\t' + str(abs(l2)) + '\n')
                    elif rand.random() <= chance:
                        rh.write('\t'.join(row.values()) + '\t' + str(l2) +
                                 '\t' + str(abs(l2)) + '\n')
"

log2="
import math
import sys
with open(sys.argv[1], 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        try:
            print(row[6], row[8], math.log(float(row[10]),2)-math.log(float(row[11]),2), sep='\t')
        except ValueError:
            continue
"

nodule_maf_filter="
import random
import sys
rand = random.Random()
l = []
with open(sys.argv[1],'rt') as handle:
    _ = handle.readline()
    for line in handle:
        maf_sum = 0.
        row = line.strip().split('\t')
        rs = row[0]
        for i in range(2, len(row)):
            af = float(row[i])
            maf_sum += min(af, 1-af)
        m = maf_sum / (len(row)-2)
        if m >= ${MIN_SYM_MAF} and 1-m >= ${MIN_SYM_MAF}:
            print(rs)
"

## Choose top unique small LD groups. Another way to do this
## would be to choose the top 200 pairs that involve a small LD
## group, but that tends to produce a shorter than desired list of genes.
extract_candidates="#!/usr/bin/env python3
import collections
import csv
import gzip
import itertools
import sys

def openfun(fname, mode):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

max_ld_groups = int(sys.argv[1])
max_ld_group_size = int(sys.argv[2])

annotation = collections.defaultdict(list)
ld_group_size = {}
betas = []
betas_small_ldgs = []
ld_groups = {}
with openfun(sys.argv[3], 'rt') as handle:
    for (ldg, beta), rows in itertools.groupby(csv.reader(handle, delimiter='\t'), lambda r: (r[1], r[3])):
        ols = set()
        for row in rows:
            if row[4] != '' and row[4] != '.': ols.update(row[4].split(','))
            annotation[ldg].append(row)
            ld_groups[row[0]] = ldg
        ld_group_size[ldg] = len(ols)
        betas.append((abs(float(row[3])), row[1], float(row[3])))

betas.sort(reverse=True)

top_medicago_hits = collections.defaultdict(set)
top_hit = {}
with open(sys.argv[4], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        afc = abs(float(row[14]))
        ldg = ld_groups[row[0]]
        if ldg in top_hit and afc > top_hit[ldg]:
            top_medicago_hits[ldg].clear()
            top_hit[ldg] = afc
            top_medicago_hits[ldg].add(row[1])
        elif ldg not in top_hit or top_hit[ldg] == afc:
            top_medicago_hits[ldg].add(row[1])
            top_hit[ldg] = afc

print('run\trs\tld_group\tp_lrt\tbeta\tortholog\tUSDA1106_ncbi\tUSDA1106_mage\tproducts\tgroup_size\tp_rank\tp_rank_small_groups\tb_rank\tb_rank_small_groups\tfdr\tbeta_with_dir\tmedicago_variants')
b_rank = 0
b_rank_sl = 0
done = set()
rs_done = set()
for b, ldg, b2 in betas:
    if ldg in done: continue
    done.add(ldg)
    b_rank += 1
    if ld_group_size[ldg] <= max_ld_group_size:
        b_rank_sl += 1
        brsl = b_rank_sl
    else:
        brsl = 'nan'
    for a in annotation[ldg]:
        if a[0] in rs_done: continue
        rs_done.add(a[0])
        print('real', a[0], a[1], 'nan', b, a[4], a[5], a[6], a[7],
              ld_group_size[ldg], 'nan', 'nan', b_rank, brsl, 'nan',
              b2, ','.join(top_medicago_hits[a[1]]), sep='\t')
"

PROJDIR="${HOME}/project/select_reseq"
GWA_PANEL="${PROJDIR}/data/strain/lists/spanfran2_medicago_truncatula.txt"  # List of accesions included in GWA panel and imputation panel
MTR_GENOTYPES="${PROJDIR}/results/combine_variants/spanfran2_medicago_2020-11-16_freebayes_rdvs/variants.bcf"
STRAINS="${PROJDIR}/data/strain/lists/spanfran2_meliloti.txt"
ENS_GENOTYPES="${PROJDIR}/results/combine_variants/meliloti_2021-04-06_denovo/variants.bcf"
PAVANNOTATION="${PROJDIR}/../jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/annotation_table.overlap.0.80.tsv.gz"
REFANNOTATION="${PROJDIR}/data/reference/USDA1106/genome.gff3"
REFANNOTATION_MAGE="${PROJDIR}/data/reference/USDA1106/genome.mage.gff3"
LD_GROUPS="${PROJDIR}/results/ld/pangenome/r2_groups/2021-04-06_denovo/SP2/${ENS_R2_THRESHOLD}/one_variant.tsv"
STRAIN_FREQUENCIES="${PROJDIR}/results/phenotypes/rhizobia_phenotypes/projected_frequency/2021-01-02/meliloti/SP2-gwas.projected_strain_frequencies.tsv"
MTR_ANNOTATION="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.gff"
MTR_ANNOTATION2="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/1.8-annotation/annotation-v1.8-ncbi-format-assembly_v1.gff"

NCR_GENES="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/ncr_search.tsv"
NCR_SCAFFOLD_NAMES="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/report.txt"

OUTDIR="${PROJDIR}/results/genome_scan/spanfran2/allele_specific_w/2021-04-07"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRATCHDIR="/scratch.global/${USER}/spanfran2_allele_specific_w_2021-04-07"

gffextract="${PROJDIR}/script/bin/extract_gff_annot.py"
w="${PROJDIR}/script/bin/intergenomic_w"
histm="${PROJDIR}/script/bin/big_histogram_multi.py"
hist_combiner="${PROJDIR}/script/bin/big_histogram_combiner.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.10.2
    module load bedtools/2.26.0
    module load htslib/1.9
    module load plink/1.90b6.10
    module load R/3.4.3

    set -euo pipefail

    datafile="${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}"
    TASK="$(cut -f 1 "$datafile")"
    BATCH="$(cut -f 2 "$datafile")"

    SECONDS=0

    cd "$SCRATCHDIR"

    case "$TASK" in

        prep)

            # Filter the Mtr SNPs
            bcftools view -Ou --samples-file "$GWA_PANEL" "$MTR_GENOTYPES" \
                | bcftools view -Ou --min-af "${MIN_HOST_MAF}:minor" \
                | bcftools view -Ov -e 'F_PASS(GT=="mis")>'"${MAX_MISS}" \
                | grep -v "rdv-" \
                > mtr_varfiltered.vcf \
                || { echo "filtering Mtr SNPs failed"; exit 1; }

            grep "^#" mtr_varfiltered.vcf > header.txt

            awk -F$'\t' '$3 == "gene" { print };' "$MTR_ANNOTATION" \
                | bedtools window -w "$GENE_WINDOW_FILTER" -a mtr_varfiltered.vcf -b - -u \
                | cat header.txt - \
                | bcftools view -Ou \
                > mtr_vargenefiltered.bcf \
                || { echo "filtering to just genic regions failed"; exit 1; }

            plink --bcf mtr_vargenefiltered.bcf \
                --allow-no-sex --allow-extra-chr \
                --indep-pairwise $THINNING_WINDOW $THINNING_STEP $THINNING_THRESHOLD \
                --out thinned \
                || { echo "making list of variants to thin from Mtr failed"; exit 1; }

            bcftools filter -Ou -i 'ID=@'"thinned.prune.in" mtr_vargenefiltered.bcf \
                > mtr_filtered.bcf \
                || { echo "thinning Mtr by LD failed"; exit 1; }

            # Filter the Ensifer SNPs and PAVs
            cut -f 4 "$LD_GROUPS" > one_per_ld_group.txt \
                || { echo "making list of representative variants failed"; exit 1; }
            bcftools view -Ov --samples-file "$STRAINS" "$ENS_GENOTYPES" \
                | bcftools filter -Ov -i 'ID=@'"one_per_ld_group.txt" \
                | bcftools view -Ov -e 'F_PASS(GT=="mis")>'"${MAX_MISS}" \
                > ens_filtered.vcf \
                || { echo "filtering Ensifer SNPs failed"; exit 1; }
            ;;

        allele-frequency)

            # Reformat strain frequency table
            sed 's/SP2_//g' "$STRAIN_FREQUENCIES" > final_freqs.tsv \
                || { echo "reformatting strain frequency header failed"; exit 1; }

            # Make table of initial frequencies
            Rscript -e 'x=read.csv("'final_freqs.tsv'", sep="\t"); x[1:nrow(x), -1] = 1/nrow(x); write.table(x, stdout(), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE);' \
                > initial_freqs.tsv \
                || { echo "reformatting strain frequency header failed"; exit 1; }

            # Calculate allele frequencies
            "$w" initial_freqs.tsv final_freqs.tsv ens_filtered.vcf \
                > ens_af.tsv \
                || { echo "calculating rhizobial allele frequencies failed"; exit 1; }

            awk -F$'\t' 'NR == 1 || $2 == "initial" { print };' ens_af.tsv \
                > ens_af.initial.tsv \
                || { echo "separating initial AF failed"; exit 1; }

            awk -F$'\t' 'NR == 1 || $2 == "final" { print };' ens_af.tsv \
                > ens_af.final.tsv \
                || { echo "separating final AF failed"; exit 1; }

            python -c "$nodule_maf_filter" ens_af.final.tsv \
                > ens_var_list.txt \
                || { echo "filtering out sites with < min MAF across nodules failed"; exit 1; }
            ;;

        run)

            split -n r/$BATCH/$N_BATCHES ens_var_list.txt > "vars.${BATCH}.txt"

            # Filter rhizobial variants more based on the amount of
            # variance in allele frequency across hosts.
            bcftools filter -Ov -i 'ID=@'"vars.${BATCH}.txt" ens_filtered.vcf \
                > "ens_filtered.${BATCH}.vcf" \
                || { echo "filtering variants for batch ${BATCH} failed"; exit 1; }

            # Run the intergenomic LD program
            "$w" initial_freqs.tsv final_freqs.tsv  "ens_filtered.${BATCH}.vcf" \
               mtr_filtered.bcf \
                > "intergenomic_ld.${BATCH}.tsv" \
                || { echo "calculating intergenomic LD failed"; exit 1; }

            # Histogram of results
            "$histm"  -n 7,9,13 0.05 0.05 0.01 "intergenomic_ld.${BATCH}.tsv" \
                > "histogram.${BATCH}.tsv" \
                || { echo "making histogram for ${BATCH} failed"; exit 1; }

            python -c "$log2" "intergenomic_ld.${BATCH}.tsv" \
                > "log2.${BATCH}.tsv" \
                || { echo "log transforming ${BATCH} failed"; exit 1; }
            "$histm"  -n 1,2,3 0.05 0.05 0.01 "log2.${BATCH}.tsv" \
                > "histogram.log2.${BATCH}.tsv" \
                || { echo "making log-transformed histogram for ${BATCH} failed"; exit 1; }
            ;;

        combine)
            "$hist_combiner" $(ls histogram.*.tsv | grep -v log2) > combined_histogram.tsv \
                || { echo "combining histograms failed"; exit 1; }
            "$hist_combiner" histogram.log2.*.tsv > combined_histogram.log2.tsv \
                || { echo "combining histograms failed"; exit 1; }
            cp combined_histogram.* "$OUTDIR"
            cp ens_var_list.txt "$OUTDIR"
            cp thinned.prune.in "$OUTDIR"
            cp final_freqs.tsv initial_freqs.tsv ens_af.tsv "$OUTDIR"
            ;;

        filter)
            python -c "$filter_candidates" "candidates.${BATCH}" \
                "intergenomic_ld.${BATCH}.tsv" \
                || { echo "extracting candidates failed"; exit 1; }
            ;;

        annotate)

#            cat candidates.*.top.tsv \
#                | sort -t$'\t' -k 15rg,15 \
#                > top_candidates.tsv
#            cat candidates.*.random.tsv \
#                | sort -t$'\t' -k 15rg,15 \
#                > random_pairs.tsv

            cd "$OUTDIR"

            ## For medicago, annotate genes
            cut -f 2,15 top_candidates.tsv \
                | sed 's/-/\t/g' \
                | awk -F$'\t' 'NF == 7 {print $4 "\t" $5-1 "\t" $5 "\t" NR ":" $1 "-" $2 "-" $3 "-" $4 "-" $5 "-" $6 ":" $7};
                     NF == 6 {print $3 "\t" $4-1 "\t" $4 "\t" NR ":" $1 "-" $2 "-" $3 "-" $4 "-" $5 ":" $6};' \
                | bedtools intersect -loj -a - -b "$MTR_ANNOTATION2" -wb -nonamecheck \
               | "$gffextract" 13 "locus_tag,product" \
                | python -c "$format_medicago" \
                > medicago_annotation.top.tsv \
                || { echo "annotating top Medicago variants failed"; exit 1; }

            ## For ensifer, add LD group and genes
            python -c "$add_ld_groups" top_candidates.tsv "$LD_GROUPS" \
                > top_candidates.ld.tsv \
                || { echo "adding LD groups to Ensifer candidates failed"; exit 1; }

            ## List all LD groups in the analysis
            python -c "$list_ld_groups" ens_var_list.txt "$LD_GROUPS" \
                > ld_groups.txt \
                || { echo "adding LD groups to Ensifer candidates failed"; exit 1; }
            
            ## Ensifer annotation
            echo "$annotate" > annotate.py
            chmod u+x annotate.py
            ./annotate.py top_candidates.ld.tsv "$PAVANNOTATION" \
                > ensifer_annotation.top.tsv

            ## Put a table of Ensifer-Medicago annotations together
            echo "$combine_partners" > combine_partners.py
            chmod u+x combine_partners.py
            ./combine_partners.py medicago_annotation.top.tsv \
                ensifer_annotation.top.tsv \
                > combined_annotation.top.tsv \
                || { echo "combining annotations failed"; exit 1; }

            ## List Ensifer candidates in the same form as the candidates
            ## from the GEMMA analyses.
            n_candidates=10000
            echo "$extract_candidates" > extract_candidates.py
            chmod u+x extract_candidates.py
            ./extract_candidates.py $n_candidates $LD_GROUP_SIZE \
                ensifer_annotation.top.tsv top_candidates.tsv \
               > ensifer_candidates.tsv \
                || { echo "making candidates file failed"; exit 1; }


            #cp final_freqs.tsv combined_annotation.top.tsv ld_groups.txt ensifer_candidates.tsv ensifer_annotation.top.tsv medicago_annotation.top.tsv top_candidates*.tsv random_pairs.tsv "$OUTDIR"
            ;;

    esac

    echo "RUN TIME = ${SECONDS} seconds or $(($SECONDS/60)) minutes"

    rm -f "$datafile"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYJOBDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-${1}-$(basename $0)"
    cp $0 $sfile

    i=0
    case "$1" in

        prep)
            WALLTIME=00:30:00
            NTHREADS=1
            MEM=64G
            echo "$1" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

        allele-frequency)
            WALLTIME=00:10:00
            NTHREADS=1
            MEM=64G
            echo "$1" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

        run)
            WALLTIME=02:00:00
            NTHREADS=1
            MEM=32G
            for b in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'$b > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        combine)
            WALLTIME=00:20:00
            NTHREADS=1
            MEM=4G
            echo "$1" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

        filter)
            WALLTIME=00:15:00
            NTHREADS=1
            MEM=4G
            for b in $(seq 1 $N_BATCHES); do
                echo "$1"$'\t'$b > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        annotate)
            WALLTIME=00:05:00
            NTHREADS=1
            MEM=8G
            echo "$1" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

    esac

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi
    array="${array}%100"

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --ntasks=$NTHREADS \
        --time=$WALLTIME --mem=$MEM --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
