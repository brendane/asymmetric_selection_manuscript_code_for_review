#!/bin/bash
#
# Run GEMMA mixed-model and linear model association analyses on SpanFran2
# fitness phenotypes, with each replicon in a separate run.
#
# Variants from denovo assemblies, one per LD group.
#
# When genes were assigned to replicons, plasmid and pSymA genes were both
# put in the pSymA bin, because it seems that genes move back and forth
# between those compartments frequently.
#
# Steps:
#   1. prep:        Create K-matrix
#   2. gwas:        Run GEMMA
#   3. gwas-lm:     Run GEMMA in plain linear model mode
#
# Rhizobia phenotype input data pipeline (a little messy):
#
# 1 relative_fitness_spanfran2_2020-10-27.r: calculate the relative fitness values for the SpanFran2 experiment
# 2 relative_fitness_compile_2020-05-23.r: put all the relative fitness values in one spot (no additional processing)
# 3 phenotypes_gwas_phenotypes_spanfran2_2020-11-13_pca.r: run PCA on relative fitness values to get community comp. phenotypes
# 4 rhizobia_phenotypes_spanfran2_2020-11-03.r: compile the median fitness values & PCA into one table, add overall average fitness
#
# SETTINGS
#
# slurm
QUEUE="amdsmall"
NTHREADS=1
NTHREADS_RANDOM=20
#
# GWAS
R2_THRESHOLD=0.95
K_MAT_TYPE=2      # 1 = centered, 2 = standardized
LMM_MODE=4        # Linear mixed model run mode, 4 = all tests (necessary for effect sizes)
MAX_MISSING=0.2   # Maximum missingness
N_PERMUTATIONS=1000
LD_GROUP_SIZE=3

TRAITS=( rhizobia_score_PC1 rhizobia_score_PC2 )
REPLICONS=( chromosome psymb psyma unclassified )

RUN="2021-06-29_irnn_by_replicon"
PROJDIR="${HOME}/project/select_reseq"
PHENOFILE="${PROJDIR}/results/phenotypes/gwas_phenotypes/spanfran2/2020-11-13_pca/strain_pca_scores.tsv"
LD_DIR="${PROJDIR}/results/ld/pangenome/r2_groups/2021-04-06_denovo"
REPLICON_FILE="${PROJDIR}/results/summary_analyses_and_data_collections/analyses/spanfran2_gwas/2021-07-30_spreadsheet/ld_group_replicons.tsv"
PAVANNOTATION="${PROJDIR}/../jgi_sequencing/results/ortho/medmel/sibeliaz/2020-04-09_split/meliloti/annotation_table.overlap.0.80.filtered.tsv.gz"

OUTDIR="${PROJDIR}/results/genome_scan/spanfran2/gwas/${RUN}"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"

awkcols2='NR==1{for(i=1;i<=NF;i++){if($i==c1){n=i}; if($i==c2){nn=i}}}; NR>1 {print $n "\t" $n "\t" $nn};'
tsv2bimbam="${PROJDIR}/script/bin/tsv2bimbam.py"

## Transformation
## Code adapted from https://yuxuanstat.com/posts/2020/06/rank-based-inverse-normal-transformation/
irrn=''
irrn+='x = scan("stdin"); '
irrn+='set.seed(61418); ' # Seed chosen randomly, just for breaking ties
irrn+='y = qnorm((rank(x, na.last="keep", ties.method="random")-0.5)/length(x)); '
irrn+='cat(y, sep="\n"); '

filter_replicons="
import csv
import sys
keep = set()
with open('${REPLICON_FILE}', 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        if row['assignment'] == sys.argv[1]:
            keep.add(row['ld_group'])
with open(sys.argv[2], 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        if 'group-' + row[4] in keep:
            sys.stdout.write(line)
"

filter_ld_groups="
import csv
import re
import sys
max_size = int(sys.argv[2])
with open(sys.argv[1], 'rt') as handle:
    for line in handle:
        row = line.strip().split('\t')
        genes = set(map(lambda v: re.sub('-.+', '', re.sub('.+-OL', 'OL', v)), row[5].split(',')))
        if len(genes) > max_size: continue
        print('\t'.join(row))
"

annotate="#!/usr/bin/env python3
import csv
import gzip
import itertools
import re
import sys
ologs = set()
annotation_info = []
with open(sys.argv[2],'rt') as handle:
    rdr = csv.reader(handle, delimiter='\t')
    for row in rdr:
        olog = re.sub('-.+', '', re.sub('.+-OL', 'OL', row[0]))
        ologs.add(olog)
        annotation_info.append([olog] + row)
with gzip.open(sys.argv[1], 'rt') as handle:
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
## variant ID, group, p-value, beta, OL, NCBI lt, MaGe lt, descriptions
for x in annotation_info:
    print(x[1], x[2], x[3], x[4], x[0], x[6], x[7], x[8], sep='\t')
"

extract_candidates="#!/usr/bin/env python3
import collections
import csv
import gzip
import itertools
import sys

N_STRAINS = 88

def openfun(fname, mode):
    if fname.endswith('.gz'):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

max_ld_groups = int(sys.argv[1])
max_ld_group_size = int(sys.argv[2])

fdr_values = {}
with openfun(sys.argv[3], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        fdr_values[row[1]] = float(row[4])

fdr_values_small = {}
with openfun(sys.argv[4], 'rt') as handle:
    for row in csv.reader(handle, delimiter='\t'):
        fdr_values_small[row[1]] = float(row[4])

annotation = collections.defaultdict(list)
ld_group_size = {}
with openfun(sys.argv[5], 'rt') as handle:
    for ldg, rows in itertools.groupby(csv.reader(handle, delimiter='\t'), lambda r: r[1]):
        ols = set()
        for row in rows:
            if row[4] != '' and row[4] != '.': ols.update(row[4].split(','))
            annotation[ldg].append(row)
        ld_group_size[ldg] = len(ols)

def af2maf(x):
    if x < 0.5:
        return x
    else:
        return 1-x

real_pves = []
real_pves_small_groups = []
with openfun(sys.argv[7], 'rt') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        eff = float(row['beta'])
        maf = af2maf(float(row['af']))
        se = float(row['se'])
        try:
            n = N_STRAINS - int(row['n_miss'])
        except KeyError:
            n = N_STRAINS - int(row['n_mis'])
        numerator = (2 * (eff**2) * maf * (1-maf))
        pve = numerator / (numerator + (se**2) * 2 * n * maf * (1-maf))
        real_pves.append((pve, row['rs']))
        if ld_group_size[row['rs']] <= max_ld_group_size:
            real_pves_small_groups.append((pve, row['rs']))
real_pves.sort(reverse=True)
real_pves_small_groups.sort(reverse=True)

print('run\trs\tld_group\tp_lrt\tbeta\tortholog\tUSDA1106_ncbi\tUSDA1106_mage\tproducts\tgroup_size\tp_rank\tp_rank_small_groups\tb_rank\tb_rank_small_groups\tfdr\tfdr_small_groups\tbeta_with_dir\tpve_rank\tpve_rank_small_groups')
with openfun(sys.argv[6], 'rt') as handle:
    for run, rrows in itertools.groupby(csv.DictReader(handle, delimiter='\t'), lambda r: r['run']):
        betas = []
        betas_small_ldgs = []
        pvals = []
        pvals_small_ldgs = []
        for row in rrows:
            betas.append((abs(float(row['beta'])), row['rs'], float(row['beta'])))
            pvals.append((float(row['p_lrt']), row['rs']))
            if ld_group_size[row['rs']] <= max_ld_group_size:
                betas_small_ldgs.append((abs(float(row['beta'])), row['rs']))
                pvals_small_ldgs.append((float(row['p_lrt']), row['rs']))
        betas.sort(reverse=True)
        pvals.sort()
        betas_small_ldgs.sort(reverse=True)
        pvals_small_ldgs.sort()
        b_ranks_sl = {ldg:i for i, (b, ldg) in enumerate(betas_small_ldgs)}
        p_ranks = {ldg:i for i, (p, ldg) in enumerate(pvals)}
        p_ranks_sl = {ldg:i for i, (p, ldg) in enumerate(pvals_small_ldgs)}
        if run == 'real':
            pve_ranks = {ldg:i for i, (pve, ldg) in enumerate(real_pves)}
            pve_ranks_sl = {ldg:i for i, (pve, ldg) in enumerate(real_pves_small_groups)}
        else:
            pve_ranks = {}; pve_ranks_sl = {}
        for b_rank, (b, ldg, b2) in enumerate(betas):
            ldg_size = ld_group_size[ldg]
            try:
                b_rank_sl = b_ranks_sl[ldg] + 1
                p_rank_sl = p_ranks_sl[ldg] + 1
                if run == 'real':
                    pve_rank_sl = pve_ranks_sl[ldg] + 1
                else:
                    pve_rank_sl = float('nan')
            except KeyError:
                b_rank_sl = float('nan')
                p_rank_sl = float('nan')
                pve_rank_sl = float('nan')
            try:
                pve_rank = pve_ranks[ldg] + 1
            except:
                pve_rank = float('nan')
            p_rank = p_ranks[ldg] + 1
            if b_rank > max_ld_groups and p_rank > max_ld_groups and b_rank_sl > max_ld_groups and p_rank_sl > max_ld_groups and run != 'real':
                continue
            if run == 'real':
                fdr = fdr_values[ldg]
                try:
                    sfdr = fdr_values_small[ldg]
                except KeyError:
                    sfdr = 'nan'
            else:
                fdr = 'nan'
                try:
                    sfdr = fdr_values[ldg]
                except KeyError:
                    sfdr = 'nan'
            for a in annotation[ldg]:
                print(run, a[0], a[1], pvals[p_rank-1][0], betas[b_rank][0],
                      '\t'.join(a[4:]), ldg_size, p_rank, p_rank_sl,
                      b_rank + 1, b_rank_sl, fdr, sfdr, b2,
                      pve_rank, pve_rank_sl, sep='\t')
"

# See Abney et al. 2015, Genetics; mvnpermute is a way of permuting phenotypes that maintains correlations
# Code for setting up covariance matrix borrowed from
# https://github.com/voichek/kmersGWAS/blob/master/src/R/transform_and_permute_phenotypes.R
# which comes from Yoav Voichek and Detlef Weigel (2019).
# Starts with the K-matrix, and then all elements are scaled by
# the genetic variance (vg). Diagonal elements also have the
# environmental variance added.
mvnpermute_code=''
mvnpermute_code+='library(mvnpermute); '
mvnpermute_code+='argv = commandArgs(trailingOnly=TRUE); '
mvnpermute_code+='K = as.matrix(read.csv(argv[1], sep="\t", header=FALSE)); '
mvnpermute_code+='phenotype = scan(argv[2]); '
mvnpermute_code+='n_permute = as.numeric(argv[3]); '
mvnpermute_code+='vge_string = scan(argv[4], what="character", sep="\n"); '
mvnpermute_code+='vg = as.numeric(strsplit(vge_string[grepl("## vg estimate", vge_string)], " ")[[1]][9]); '
mvnpermute_code+='ve = as.numeric(strsplit(vge_string[grepl("## ve estimate", vge_string)], " ")[[1]][9]); '
mvnpermute_code+='covariance_matrix = vg * K + ve * diag(dim(K)[1]); '
mvnpermute_code+='permute_phenotype = mvnpermute(phenotype, rep(1, length(phenotype)), covariance_matrix, nr=n_permute); '
mvnpermute_code+='write.table(permute_phenotype, file=stdout(), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '

# Completely randomized phenotypes
permute_code=''
permute_code+='argv = commandArgs(trailingOnly=TRUE); '
permute_code+='phenotype = scan(argv[1]); '
permute_code+='n_permute = as.numeric(argv[2]); '
permute_code+='permute_phenotype = matrix(nrow=length(phenotype), ncol=n_permute); '
permute_code+='for(i in 1:n_permute) { '
permute_code+='permute_phenotype[, i] = sample(phenotype, length(phenotype), replace=FALSE); }; '
permute_code+='write.table(permute_phenotype, file=stdout(), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE); '


if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.9
    module load gemma/0.98.1
    module load parallel/20190122
    module load R/4.0.4

    set -euo pipefail

    datafile="${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}"
    TASK="$(cut -f 1 "$datafile")"
    PHENOTYPE="$(cut -f 2 "$datafile")"

    SECONDS=0

    cd "$OUTDIR"

    variants="${LD_DIR}/SP2/${R2_THRESHOLD}/filtered_variants.tsv.gz"
    onevar="${LD_DIR}/SP2/${R2_THRESHOLD}/one_variant.tsv"

    case "$TASK" in

        prep)
            # Copy strains from first column of phenotype file to make
            # sure that the list of strains will be in the same order
            # as the phenotypes
            cut -f 1 "$PHENOFILE" | tail -n +2 | sort > strains.txt \
                || { echo "copying strains file failed"; exit 1; }

            # Fake phenotype file to make gemma happy
            rm -f variants.p
            while read -r strain; do
                echo 1 >> variants.p
            done < strains.txt

            for replicon in "${REPLICONS[@]}"; do
                # Make a list of LD groups that are part of this replicon.
                # Use the majority-rule analysis that was used for the candidate
                # gene comparisons.
                python3 -c "$filter_replicons" "$replicon" "$onevar" > \
                    "onevar.${replicon}" \
                    || { echo "filtering by replicon (${replicon}) failed"; exit 1; }

                # Convert to bimbam format. Do this in a way that ensures
                # the genotypes are in the same order as the phenotype files
                # by taking the order of strains from the phenotype table and
                # using --order-by-keep.
                cp "$tsv2bimbam" .
                "./$(basename "$tsv2bimbam")" \
                    --use-2 \
                    --one-variant "onevar.${replicon}" \
                    --keep strains.txt \
                    --order-by-keep \
                    "$variants" \
                    | sed 's/nan/NA/g' \
                    > "variants.${replicon}.g" \
                    || { echo "converting to bimbam format for failed"; exit 1; }

                # Calculate relatedness matrix
                gemma -maf 0.01 -miss "$MAX_MISSING" \
                    -g "variants.${replicon}.g" \
                    -p variants.p \
                    -gk "$K_MAT_TYPE" \
                    -o "K.${replicon}" \
                    || { echo "calculating K matrix for failed"; exit 1; }
            done

            # Delete intermediate files for K-matrix calculation
            rm -rf "K"
            mv "output" "K"
            ;;


        mvnpermute)
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            # Create a phenotype input file for gemma from real data
            # Manually check that order of phenotype values is correct;
            # Should match the order of the genotypes file and strains.txt.
            rm -rf phenotypes
            mkdir -p phenotypes
            awk -F$'\t' -v c1="strain" -v c2="$PHENOTYPE" "$awkcols2" \
                "$PHENOFILE" \
                | sort -t$'\t' -k 1b,1 \
                | join -t$'\t' -j 1 ../strains.txt - \
                | cut -f 3 \
                | Rscript -e "$irrn" \
                > phenotypes/real.p \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            for replicon in "${REPLICONS[@]}"; do
                ## -lmm 5: estimate Vg and Ve
                gemma -lmm 5 \
                    -maf 0.01 \
                    -miss "$MAX_MISSING" \
                    -o "vge.${replicon}" \
                    -g "../variants.${replicon}.g" \
                    -p phenotypes/real.p \
                    -k "../K/K.${replicon}.sXX.txt" \
                    || { echo "estimating variance components for ${PHENOTYPE} failed"; exit 1; }

                # Run MVNPermute on this phenotype
                Rscript -e "$mvnpermute_code" "../K/K.${replicon}.sXX.txt" \
                    phenotypes/real.p $N_PERMUTATIONS \
                    "output/vge.${replicon}.log.txt" \
                    > phenotypes/mvn_random.${replicon}.tsv \
                    || { echo "creating random phenotypes ${PHENOTYPE} failed"; exit 1; }
            done
            ;;



        gwas)
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            for replicon in "${REPLICONS[@]}"; do

                variants="../variants.${replicon}.g"
                kmat="../K/K.${replicon}.sXX.txt"
                comb_output="combined_output.${replicon}.tsv"
                comb_log="combined_log.${replicon}.txt"
                real_prefix="output.real.${replicon}"

                gemma \
                    -maf 0.01 \
                    -miss "$MAX_MISSING" \
                    -lmm "$LMM_MODE" \
                    -g "$variants" \
                    -p phenotypes/real.p \
                    -k "$kmat" \
                    -o "$real_prefix" \
                    || { echo "gemma -lmm 4 failed on real data, ${PHENOTYPE}"; exit 1; }

                awk -F$'\t' 'NR == 1 {print "run\t" $2 "\t" $7 "\t" $8 "\t" $9 "\t" $14}; NR > 1 {print "real\t" $2 "\t" $7 "\t" $8 "\t" $9 "\t" $14};' \
                    "output/${real_prefix}.assoc.txt" \
                    | sort -t$'\t' -k 6g,6 \
                    > "$comb_output" \
                    || { echo "print real gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }

                echo "### real" > "$comb_log"
                cat "output/${real_prefix}.log.txt" >> "$comb_log" \
                    || { echo "print real gemma log for ${PHENOTYPE} to combined file failed"; exit 1; }
                gzip -f "output/${real_prefix}.assoc.txt" \
                    || { echo "gzipping gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }
                rm "output/${real_prefix}.log.txt"

                rungemma() {
                    iter=$1
                    # Create a phenotype input file for plink.
                    cut -f $(($iter+1)) "phenotypes/mvn_random.${replicon}.tsv" \
                        > "phenotype.${iter}.p" \
                        || { echo "making random phenotype file for ${iter}, ${PHENOTYPE} failed"; exit 1; }

                    # Run GEMMA's simplest LMM
                    # -lmm 4 = all p-value tests - necessary to get the beta estimate
                    # -k = Use the k matrix made in the first step
                    gemma \
                        -maf 0.01 \
                        -miss "$MAX_MISSING" \
                        -g "$variants" \
                        -p "phenotype.${iter}.p" \
                        -k "$kmat" \
                        -lmm "$LMM_MODE" \
                        -o "output.random.${iter}" \
                        || { echo "gemma failed on iter ${iter}, ${PHENOTYPE}"; exit 1; }
                }
                export -f rungemma

                source "$(which env_parallel.bash)"
                seq 0 $(($N_PERMUTATIONS-1)) \
                    | env_parallel -j$NTHREADS_RANDOM rungemma \
                    || { echo "random run ${PHENOTYPE} failed"; exit 1; }

                rm -f top_ld_groups.tsv
                for iter in $(seq 0 $(($N_PERMUTATIONS-1))); do
                    awk -F$'\t' 'NR > 1 {print "random-'$iter'\t" $2 "\t" $7 "\t" $8 "\t" $9 "\t" $14}' \
                        "output/output.random.${iter}.assoc.txt" \
                        | sort -t$'\t' -k 6g,6 \
                        | tee tmp \
                        >> "$comb_output" \
                        || { echo "print random gemma output for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                    awk -F$'\t' 'NR == 1 { print $1 "\t" $6 };' tmp >> top_ld_groups.tsv \
                        || { echo "getting top LD group p-value for ${iter}, ${PHENOTYPE} failed"; exit 1; }

                    echo "### random-${iter}" >> "$comb_log"
                    cat "output/output.random.${iter}.log.txt" >> "$comb_log" \
                        || { echo "print random gemma log for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                    rm "phenotype.${iter}.p" "output/output.random.${iter}."* tmp
                done

                gzip -f "$comb_output"

            done

            ;;


        gwas-lm)
            cd "$PHENOTYPE"
            mkdir -p lm && cd lm


            for replicon in "${REPLICONS[@]}"; do

                variants="../../variants.${replicon}.g"
                comb_output="combined_output.${replicon}.tsv"
                comb_log="combined_log.${replicon}.txt"
                real_prefix="output.real.${replicon}"

                gemma \
                    -maf 0.01 \
                    -miss "$MAX_MISSING" \
                    -lm "$LMM_MODE" \
                    -g "$variants"  \
                    -p ../phenotypes/real.p \
                    -o "$real_prefix" \
                    || { echo "gemma -lmm 4 failed on real data, ${PHENOTYPE}"; exit 1; }

                awk -F$'\t' 'NR == 1 {print "run\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12}; NR > 1 {print "real\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12};' \
                    "output/${real_prefix}.assoc.txt" \
                    | sort -t$'\t' -k 6g,6 \
                    > "$comb_output" \
                    || { echo "print real gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }

                echo "### real" > "$comb_log"
                cat "output/${real_prefix}.log.txt" >> "$comb_log" \
                    || { echo "print real gemma log for ${PHENOTYPE} to combined file failed"; exit 1; }
                gzip -f "output/${real_prefix}.assoc.txt" \
                    || { echo "gzipping gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }
                rm "output/${real_prefix}.log.txt"

                rungemma() {
                    iter=$1
                    # Create a phenotype input file for plink.
                    cut -f $(($iter+1)) "../phenotypes/mvn_random.${replicon}.tsv" \
                        > "phenotype.${iter}.p" \
                        || { echo "making random phenotype file for ${iter}, ${PHENOTYPE} failed"; exit 1; }

                    # Run GEMMA's simplest LMM
                    # -lmm 4 = all p-value tests - necessary to get the beta estimate
                    # -k = Use the k matrix made in the first step
                    gemma \
                        -maf 0.01 \
                        -miss "$MAX_MISSING" \
                        -g "$variants" \
                        -p "phenotype.${iter}.p" \
                        -k "$kmat" \
                        -lm "$LMM_MODE" \
                        -o "output.random.${iter}" \
                        || { echo "gemma failed on iter ${iter}, ${PHENOTYPE}"; exit 1; }
                }
                export -f rungemma

                source "$(which env_parallel.bash)"
                seq 0 $(($N_PERMUTATIONS-1)) \
                    | env_parallel -j$NTHREADS_RANDOM rungemma \
                    || { echo "random run ${PHENOTYPE} failed"; exit 1; }

                rm -f top_ld_groups.tsv
                for iter in $(seq 0 $(($N_PERMUTATIONS-1))); do
                    awk -F$'\t' 'NR > 1 {print "random-'$iter'\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12};' \
                        "output/output.random.${iter}.assoc.txt" \
                        | sort -t$'\t' -k 6g,6 \
                        | tee tmp \
                        >> "$comb_output" \
                        || { echo "print random gemma output for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                    awk -F$'\t' 'NR == 1 { print $1 "\t" $6 };' tmp >> top_ld_groups.tsv \
                        || { echo "getting top LD group p-value for ${iter}, ${PHENOTYPE} failed"; exit 1; }

                    echo "### random-${iter}" >> "$comb_log"
                    cat "output/output.random.${iter}.log.txt" >> "$comb_log" \
                        || { echo "print random gemma log for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                    rm "phenotype.${iter}.p" "output/output.random.${iter}."* tmp
                done

                gzip -f "$comb_output"

            done

            ;;


        annotate*)
            cd "$PHENOTYPE"
            if [[ "$TASK" == "annotate-lm" ]]; then
                cd "lm"
            fi

            for replicon in "${REPLICONS[@]}"; do
                onevar="../onevar.${replicon}"
                if [[ "$TASK" == "annotate-lm" ]]; then
                    onevar="../../onevar.${replicon}"
                fi

                comb_output="combined_output.${replicon}.tsv.gz"

                mkdir -p "annotation-${replicon}"

                awk -F$'\t' 'BEGIN {OFS="\t"}; {$5 = "group-" $5; print $0};' \
                    "$onevar" > tmp.ov \
                    || { echo "tweaking one variant file for ${PHENOTYPE} failed"; exit 1; }

                "${PROJDIR}/script/bin/gemma_significant_from_combined_output.py" \
                    --fdr 0.01,0.05,0.10,0.20,1 \
                    --effect-sizes \
                    --one-variant tmp.ov "$comb_output" \
                    > "annotation-${replicon}/all_groups.tsv" \
                    || { echo "getting list of significant variants for ${PHENOTYPE} failed"; exit 1; }

                python3 -c "$filter_ld_groups" tmp.ov $LD_GROUP_SIZE > small.ov \
                    || { echo "filtering to small LD groups for ${PHENOTYPE} failed"; exit 1; }
                "${PROJDIR}/script/bin/gemma_significant_from_combined_output.py" \
                    --only-known-groups \
                    --fdr 0.01,0.05,0.10,0.20,1 \
                    --effect-sizes \
                    --one-variant small.ov "$comb_output" \
                    > "annotation-${replicon}/all_groups.small.tsv" \
                    || { echo "getting list of significant variants in small LD groups for ${PHENOTYPE} failed"; exit 1; }
                rm tmp.ov small.ov

                ## Annotation
                echo "$annotate" > annotate.py
                chmod u+x annotate.py
                ./annotate.py "$PAVANNOTATION" \
                    "annotation-${replicon}/all_groups.tsv" \
                    > "annotation-${replicon}/all.annotation.tsv" \
                    || { echo "getting annotations for ${PHENOTYPE} failed"; exit 1; }

                n_candidates=400
                echo "$extract_candidates" > extract_candidates.py
                chmod u+x extract_candidates.py
                ./extract_candidates.py \
                    $n_candidates \
                    $LD_GROUP_SIZE \
                    "annotation-${replicon}/all_groups.tsv" \
                    "annotation-${replicon}/all_groups.small.tsv" \
                    "annotation-${replicon}/all.annotation.tsv" \
                    "$comb_output" \
                    "output/output.real.${replicon}.assoc.txt.gz" \
                    > "annotation-${replicon}/candidates.tsv" \
                    || { echo "making candidates file for ${PHENOTYPE} failed"; exit 1; }

                grep "^run\|^real" "annotation-${replicon}/candidates.tsv" \
                    > "annotation-${replicon}/real_candidates.tsv" \
                    || { echo "extracting real data for ${PHENOTYPE} failed"; exit 1; }
                gzip -f "annotation-${replicon}/candidates.tsv" \
                    || { echo "gzipping ${PHENOTYPE} failed"; exit 1; }
                gzip -f "annotation-${replicon}/real_candidates.tsv" \
                    || { echo "gzipping ${PHENOTYPE} failed"; exit 1; }
            done
            ;;

    esac

    # See how long the script really took for future reference
    echo "RUN TIME ${PHENOTYPE} = ${SECONDS} seconds or $(($SECONDS/60)) minutes"

    rm -f $datafile

else

    mkdir -p "$OUTDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYJOBDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-${1}-$(basename $0)"
    cp $0 $sfile

    # Run each phenotype in a separate job (of a single array job)
    i=0
    case $1 in

        prep)
            WALLTIME="00:15:00"
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=1
            ;;

        mvnpermute)
            WALLTIME="00:05:00"
            for pheno in ${TRAITS[@]}; do
                echo "$1"$'\t'"$pheno" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        gwas)
            WALLTIME="01:00:00"
            NTHREADS=$NTHREADS_RANDOM
            for pheno in ${TRAITS[@]}; do
                echo "$1"$'\t'"$pheno" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        gwas-lm)
            WALLTIME="01:00:00"
            NTHREADS=$NTHREADS_RANDOM
            for pheno in ${TRAITS[@]}; do
                echo "$1"$'\t'"$pheno" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        #annotate*)
        annotate)
            WALLTIME="01:30:00"
            NTHREADS=1
            for task in annotate annotate-lm; do
                for pheno in ${TRAITS[@]}; do
                    echo "$task"$'\t'"$pheno" > "${ARRAYJOBDIR}/${i}"
                    i=$(($i+1))
                done
            done
            ;;

    esac

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --ntasks-per-node=$NTHREADS \
        --time=$WALLTIME --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
