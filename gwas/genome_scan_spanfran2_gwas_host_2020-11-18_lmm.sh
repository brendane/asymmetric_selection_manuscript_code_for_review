#!/bin/bash
#
# Run GEMMA's univariate LMM on the host PCA community phenotypes
# from the winter 2020 (SpanFran2) experiment. And also added the other
# host phenotypes.
#
# UPDATE: 2021-04-06: Regular linear model run added.
# UPDATE 2021-05-27: Re-ran PC traits with updated strain fitness values,
#      removed "gwas" step b/c it is redundant with gwas-effects
#
# NOTE 2021-06-29: nodule traits run in new script (2021-06-29_nod_traits_lmm),
# unadjusted means
#
# SETTINGS
#
# pbs
QUEUE=amdsmall
#
# GWAS
K_MAT_TYPE=2      # 1 = centered, 2 = standardized
LMM_MODE=2        # Linear mixed model run mode, 4 = all tests (necessary for effect sizes), 2 = just LRT (faster than all tests)
MIN_MAF=0.05      # Filter applied for all individuals and for each phenotype
MAX_MISS=0.2
N_PERMUTATIONS=100
N_PERMUTATIONS_RUN=100
NTHREADS_RANDOM=10

TRAITS=( host_score_PC1 host_score_PC2 )
PLANT_ONLY_TRAITS=( c86_hght_wk3 c86_flow_date c86_biomass c86_nods c86_nod_mass c86_nod_area c86_nod_circ c86_nod_branch biomass_norm_diff )
## Don't run permutations on biomass_norm_diff because pve was so low
PLANT_ONLY_TRAITS=( c86_hght_wk3 c86_flow_date c86_biomass c86_nods c86_nod_mass c86_nod_area c86_nod_circ c86_nod_branch )

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


PROJDIR="${HOME}/project/select_reseq"
PHENOFILE="${PROJDIR}/results/phenotypes/mtr_phenotypes/spanfran2/2020-11-03/phenotypes.tsv"
PLANT_ONLY_PHENOFILE="${PROJDIR}/results/phenotypes/mtr_phenotypes/spanfran2/2020-08-12/medians.tsv"
GWA_PANEL="${PROJDIR}/data/strain/lists/spanfran2_medicago_truncatula.txt"  # List of accesions included in GWA panel and imputation panel
GENOTYPES="${PROJDIR}/results/combine_variants/spanfran2_medicago_2020-11-16_freebayes_rdvs/variants.bcf"
ANNOTATION="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.gff"
FAI="${PROJDIR}/../mtr_variant_calls/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.fasta.fai"

OUTDIR="${PROJDIR}/results/genome_scan/spanfran2/gwas/host/2020-11-18_lmm"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRATCHDIR="/scratch.global/${USER}/spanfran2_lmm_2020-11-18"

awkcols2='NR==1{for(i=1;i<=NF;i++){if($i==c1){n=i}; if($i==c2){nn=i}}}; NR>1 {print $n "\t" $n "\t" $nn};'
manhattan="${PROJDIR}/script/bin/manhattan_plot.r"
thresh="${PROJDIR}/script/bin/gwas_empirical_threshold.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bedtools/2.29.2
    module load gcta/1.93.2
    module load gemma/0.98.1
    module load plink/1.90b6.10
    module load bcftools/1.10.2
    module load parallel
    module unload R; module load R/3.4.3

    set -euo pipefail

    datafile="${ARRAYJOBDIR}/${SLURM_ARRAY_TASK_ID}"
    TASK="$(cut -f 1 "$datafile")"
    PHENOTYPE="$(cut -f 2 "$datafile")"
    PHENOTYPE_TABLE="$(cut -f 3 "$datafile")"
    RERUN="$(cut -f 4 "$datafile")"

    SECONDS=0

    cd "$OUTDIR"

    case "$TASK" in

        prep)
            # Need some phenotypes to stick in the plink file
            cat "$GWA_PANEL" | \
                Rscript -e 'x=scan("stdin", what="character"); write.table(cbind(x, x, rnorm(length(x))), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE);' > \
                ph.txt \
                || { echo "getting random phenotypes to make plink happy failed"; exit 1; }

            # Convert to plink format
            # For this test run, just eliminate sites that are
            # biallelic in the GWA panel. At some point, these
            # should probably be reformatted as multiple variants.
            # Filter with bcftools first, because I think this uses
            # less memory.
            bcftools view -Ou --samples-file "$GWA_PANEL" "$GENOTYPES" \
                | bcftools view -Ou --min-af "${MIN_MAF}:minor" \
                | bcftools view -Ov -e 'F_PASS(GT=="mis")>'"${MAX_MISS}" \
                > filtered.vcf

            plink --vcf filtered.vcf --make-bed --out plink --double-id \
                --allow-extra-chr --allow-no-sex --memory 64000 \
                --pheno ph.txt \
                || { echo "making plink file failed"; exit 1; }

            # Calculate relatedness matrix
            gemma -bfile plink -gk $K_MAT_TYPE -o K \
                || { echo "calculating K matrix failed"; exit 1; }

            # Move output files
            rm -rf K filtered.vcf
            mv output K
            ;;

        mvnpermute)
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"
            kmat=../K/K.sXX.txt

            # Create a phenotype input file for gemma from real data
            # Manually check that order of phenotype values is correct;
            # Should match the order of the genotypes file and strains.txt.
            rm -rf phenotypes
            mkdir -p phenotypes

            strain_col=strain
            grep '^strain' "$PHENOTYPE_TABLE" || strain_col=genotype

            awk -F$'\t' -v c1="$strain_col" -v c2="$PHENOTYPE" "$awkcols2" \
                "$PHENOTYPE_TABLE" > phenotype.tsv \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            # Check that there are no real -9s in the phenotype data (used as
            # missing code by plink).
            awk -F$'\t' '$3 == -9 {print};' phenotype.tsv | grep "." && \
                { echo "overlap with missing data code in ${PHENOTYPE}: failed"; exit 1; }

            # Make a binary plink file with this phenotype included
            plink --bfile ../plink --allow-extra-chr \
                --memory 20000 \
                --make-bed --allow-no-sex --out "plink" \
                --pheno "phenotype.tsv" \
                || { echo "making plink file for ${PHENOTYPE} failed"; exit 1; }

            awk '{ print $6 };' plink.fam > phenotype.p \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            ## -lmm 5: estimate Vg and Ve
            gemma -lmm 5 \
                -o vge \
                -p phenotype.p \
                -k "$kmat" \
                || { echo "estimating variance components for ${PHENOTYPE} failed"; exit 1; }

            # Run MVNPermute on this phenotype
            Rscript -e "$mvnpermute_code" ../K/K.sXX.txt \
                phenotype.p $N_PERMUTATIONS \
                output/vge.log.txt \
                > phenotypes/mvn_random.tsv \
                || { echo "creating random phenotypes ${PHENOTYPE} in ${COMMUNITY} failed"; exit 1; }

            rm plink*
            ;;

        gwas-effects)
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            kmat="${OUTDIR}/K/K.sXX.txt"
            real_prefix=output.real.effects

            strain_col=strain
            grep '^strain' "$PHENOTYPE_TABLE" || strain_col=genotype

            # Create a phenotype input file for plink.
            awk -F$'\t' -v c1="$strain_col" -v c2="$PHENOTYPE" "$awkcols2" \
                "$PHENOTYPE_TABLE" > phenotype.tsv \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            # Check that there are no real -9s in the phenotype data (used as
            # missing code by plink).
            awk -F$'\t' '$3 == -9 {print};' phenotype.tsv | grep "." && \
                { echo "overlap with missing data code: failed"; exit 1; }

            # Make a binary plink file with this phenotype included
            plink --bfile "${OUTDIR}/plink" --allow-extra-chr \
                --memory 5000 \
                --make-bed --allow-no-sex --out "plink" \
                --pheno "phenotype.tsv" \
                || { echo "making plink file for ${PHENOTYPE} failed"; exit 1; }

            # Run GEMMA
            gemma \
                -maf $MIN_MAF \
                -miss $MAX_MISS \
                -lmm 4 \
                -bfile plink \
                -k "$kmat" \
                -o "$real_prefix" \
                || { echo "gemma failed on real data, ${PHENOTYPE}"; exit 1; }

            rm plink*
            ;;

        gcta)
            cd "$PHENOTYPE"

            cut -f 2 ../plink.bim \
                | grep -v "rdv" \
                > snps.txt \
                || { echo "extracting SNPs for ${PHENOTYPE} failed"; exit 1; }

            # Make a binary plink file with this phenotype included
            plink --bfile "${OUTDIR}/plink" --allow-extra-chr \
                --memory 5000 \
                --make-bed --allow-no-sex --out "plink" \
                --geno $MAX_MISS \
                --extract snps.txt \
                --pheno phenotype.tsv \
                || { echo "making plink file for ${PHENOTYPE} failed"; exit 1; }

            sed 's/^CM010648.1/1/g' plink.bim \
                | sed 's/^CM010649.1/2/g' \
                | sed 's/^CM010650.1/3/g' \
                | sed 's/^CM010651.1/4/g' \
                | sed 's/^CM010652.1/5/g' \
                | sed 's/^CM010653.1/6/g' \
                | sed 's/^CM010654.1/7/g' \
                | sed 's/^CM010655.1/8/g' \
                | sed 's/^PSQE01000009.1/9/g' \
                | sed 's/^PSQE01000010.1/10/g' \
                | sed 's/^PSQE01000011.1/11/g' \
                | sed 's/^PSQE01000012.1/12/g' \
                | sed 's/^PSQE01000013.1/13/g' \
                | sed 's/^PSQE01000014.1/14/g' \
                | sed 's/^PSQE01000015.1/15/g' \
                | sed 's/^PSQE01000016.1/16/g' \
                | sed 's/^PSQE01000017.1/17/g' \
                | sed 's/^PSQE01000018.1/18/g' \
                | sed 's/^PSQE01000019.1/19/g' \
                | sed 's/^PSQE01000020.1/20/g' \
                | sed 's/^PSQE01000021.1/21/g' \
                | sed 's/^PSQE01000022.1/22/g' \
                | sed 's/^PSQE01000023.1/23/g' \
                | sed 's/^PSQE01000024.1/24/g' \
                | sed 's/^PSQE01000025.1/25/g' \
                | sed 's/^PSQE01000026.1/26/g' \
                | sed 's/^PSQE01000027.1/27/g' \
                | sed 's/^PSQE01000028.1/28/g' \
                | sed 's/^PSQE01000029.1/29/g' \
                | sed 's/^PSQE01000030.1/30/g' \
                | sed 's/^PSQE01000031.1/31/g' \
                > tmp \
                && mv tmp plink.bim \
                || { echo "fixing chromosome names for ${PHENOTYPE} failed"; exit 1; }

            gcta64 --bfile plink --autosome-num 31 --maf $MIN_MAF --make-grm \
                --out gcta_output --thread-num 1 \
                || { echo "making GRM for ${PHENOTYPE} failed"; exit 1; }

            gcta64 --grm gcta_output --pheno phenotype.tsv --reml \
                --out gcta_output --thread-num 1 \
                || { echo "calculating heritability for ${PHENOTYPE} failed"; exit 1; }

            #rm plink*
            ;;

        gwas-permuted)
            cd "$SCRATCHDIR"
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            variants="${OUTDIR}/plink"
            kmat="${OUTDIR}/K/K.sXX.txt"
            comb_output=combined_output.tsv
            comb_log=combined_log.txt
            real_prefix=output.real.effects

            awk -F$'\t' 'NR == 1 {print "run\t" $2 "\t" $8 "\t" $14}; NR > 1 {print "real\t" $2 "\t" $8 "\t" $14};' \
                "${OUTDIR}/${PHENOTYPE}/output/${real_prefix}.assoc.txt" \
                | sort -t$'\t' -k 3g,3 \
                > "$comb_output" \
                || { echo "print real gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }

            echo "### real" > "$comb_log"
            cat "${OUTDIR}/${PHENOTYPE}/output/${real_prefix}.log.txt" >> "$comb_log" \
                || { echo "print real gemma log for ${PHENOTYPE} to combined file failed"; exit 1; }

            rungemma() {
                iter=$1

                if [[ "$RERUN" == "no" ]]; then
                    if [[ -f "output/output.random.${iter}.assoc.txt" && ! -f "plink.${iter}.bed" ]]; then
                        echo "skipping ${iter} because it has already been run"
                        return 0
                    fi
                fi

                # Create a phenotype input file and convert to plink format
                cut -f $(($iter+1)) "${OUTDIR}/${PHENOTYPE}/phenotypes/mvn_random.tsv" \
                    | paste "${OUTDIR}/plink.fam" - \
                    | awk '{ print $1 "\t" $1 "\t" $7 };' \
                    > "phenotype.${iter}.p" \
                    || { echo "making random phenotype file for ${iter}, ${PHENOTYPE} failed"; exit 1; }

                plink --bfile "${OUTDIR}/plink" --allow-extra-chr \
                    --memory 5000 \
                    --make-bed --allow-no-sex --out "plink.${iter}" \
                    --pheno "phenotype.${iter}.p" \
                    || { echo "making plink file for ${PHENOTYPE}, ${iter} failed"; exit 1; }

                # Run GEMMA's simplest LMM
                # -k = Use the k matrix made in the first step
                gemma \
                    -maf $MIN_MAF \
                    -miss $MAX_MISS \
                    -bfile "plink.${iter}" \
                    -k "$kmat" \
                    -lmm 4 \
                    -o "output.random.${iter}" \
                    || { echo "gemma failed on iter ${iter}, ${PHENOTYPE}"; exit 1; }

                rm "plink.${iter}."* "phenotype.${iter}.p"
            }
            export -f rungemma

            source "$(which env_parallel.bash)"
            seq 0 $(($N_PERMUTATIONS_RUN-1)) \
                | env_parallel -j$NTHREADS_RANDOM rungemma \
                || { echo "random run ${PHENOTYPE} failed"; exit 1; }

            for iter in $(seq 0 $(($N_PERMUTATIONS_RUN-1))); do
                awk -F$'\t' 'NR > 1 {print "random-'$iter'\t" $2 "\t" $8 "\t" $14};' \
                    "output/output.random.${iter}.assoc.txt" \
                    | sort -t$'\t' -k 3g,3 \
                    | tee tmp \
                    >> "$comb_output" \
                    || { echo "print random gemma output for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                echo "### random-${iter}" >> "$comb_log"
                cat "output/output.random.${iter}.log.txt" >> "$comb_log" \
                    || { echo "print random gemma log for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }

                #rm "output/output.random.${iter}."* tmp
            done

            gzip -f "$comb_output"
            ;;


        lm)
            mkdir -p lm && cd lm
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            kmat="${OUTDIR}/K/K.sXX.txt"
            real_prefix=output.real.effects

            strain_col=strain
            grep '^strain' "$PHENOTYPE_TABLE" || strain_col=genotype

            # Create a phenotype input file for plink.
            awk -F$'\t' -v c1="$strain_col" -v c2="$PHENOTYPE" "$awkcols2" \
                "$PHENOTYPE_TABLE" > phenotype.tsv \
                || { echo "making phenotype file for ${PHENOTYPE} failed"; exit 1; }

            # Check that there are no real -9s in the phenotype data (used as
            # missing code by plink).
            awk -F$'\t' '$3 == -9 {print};' phenotype.tsv | grep "." && \
                { echo "overlap with missing data code: failed"; exit 1; }

            # Make a binary plink file with this phenotype included
            plink --bfile "${OUTDIR}/plink" --allow-extra-chr \
                --memory 5000 \
                --make-bed --allow-no-sex --out "plink" \
                --pheno "phenotype.tsv" \
                || { echo "making plink file for ${PHENOTYPE} failed"; exit 1; }

            # Run GEMMA
            gemma \
                -maf $MIN_MAF \
                -miss $MAX_MISS \
                -lm 4 \
                -bfile plink \
                -o "$real_prefix" \
                || { echo "gemma failed on real data, ${PHENOTYPE}"; exit 1; }

            rm plink*
            ;;


        lm-permuted)
            cd "$SCRATCHDIR"
            mkdir -p lm && cd lm
            mkdir -p "$PHENOTYPE" && cd "$PHENOTYPE"

            variants="${OUTDIR}/plink"
            kmat="${OUTDIR}/K/K.sXX.txt"
            comb_output=combined_output.tsv
            comb_log=combined_log.txt
            real_prefix=output.real.effects

            ## TODO: Check column numbers
            awk -F$'\t' 'NR == 1 {print "run\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12}; NR > 1 {print "real\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12};' \
                "${OUTDIR}/lm/${PHENOTYPE}/output/${real_prefix}.assoc.txt" \
                | sort -t$'\t' -k 3g,3 \
                > "$comb_output" \
                || { echo "print real gemma output for ${PHENOTYPE} to combined file failed"; exit 1; }

            echo "### real" > "$comb_log"
            cat "${OUTDIR}/lm/${PHENOTYPE}/output/${real_prefix}.log.txt" >> "$comb_log" \
                || { echo "print real gemma log for ${PHENOTYPE} to combined file failed"; exit 1; }

#            rungemma() {
#                iter=$1
#
#                if [[ "$RERUN" == "no" ]]; then
#                    if [[ -f "output/output.random.${iter}.assoc.txt" && ! -f "plink.${iter}.bed" ]]; then
#                        return
#                    fi
#                fi
#
#                # Create a phenotype input file and convert to plink format
#                cut -f $(($iter+1)) "${OUTDIR}/${PHENOTYPE}/phenotypes/mvn_random.tsv" \
#                    | paste "${OUTDIR}/plink.fam" - \
#                    | awk '{ print $1 "\t" $1 "\t" $7 };' \
#                    > "phenotype.${iter}.p" \
#                    || { echo "making random phenotype file for ${iter}, ${PHENOTYPE} failed"; exit 1; }
#
#                plink --bfile "${OUTDIR}/plink" --allow-extra-chr \
#                    --memory 5000 \
#                    --make-bed --allow-no-sex --out "plink.${iter}" \
#                    --pheno "phenotype.${iter}.p" \
#                    || { echo "making plink file for ${PHENOTYPE}, ${iter} failed"; exit 1; }
#
#                # Run GEMMA's simplest LMM
#                gemma \
#                    -maf $MIN_MAF \
#                    -miss $MAX_MISS \
#                    -bfile "plink.${iter}" \
#                    -lm 4 \
#                    -o "output.random.${iter}" \
#                    || { echo "gemma failed on iter ${iter}, ${PHENOTYPE}"; exit 1; }
#
#                rm "plink.${iter}."* "phenotype.${iter}.p"
#            }
#            export -f rungemma
#
#            source "$(which env_parallel.bash)"
#            seq 0 $(($N_PERMUTATIONS_RUN-1)) \
#                | env_parallel -j$NTHREADS_RANDOM rungemma \
#                || { echo "random run ${PHENOTYPE} failed"; exit 1; }
#
#            for iter in $(seq 0 $(($N_PERMUTATIONS_RUN-1))); do
#                awk -F$'\t' 'NR > 1 {print "random-'$iter'\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $12};' \
#                    "output/output.random.${iter}.assoc.txt" \
#                    | sort -t$'\t' -k 3g,3 \
#                    | tee tmp \
#                    >> "$comb_output" \
#                    || { echo "print random gemma output for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }
#
#                echo "### random-${iter}" >> "$comb_log"
#                cat "output/output.random.${iter}.log.txt" >> "$comb_log" \
#                    || { echo "print random gemma log for ${iter}, ${PHENOTYPE} to combined file failed"; exit 1; }
#
#                rm "output/output.random.${iter}."* tmp
#            done

            gzip -f "$comb_output"
            ;;


        annotate2)
            comb_output=combined_output.tsv.gz

            ## Annotation for full set of results
            cd "$SCRATCHDIR/${PHENOTYPE}"
            cp "$comb_output" "${OUTDIR}/${PHENOTYPE}"

#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 14g,14 -t$'\t' \
#                | awk -F$'\t' 'NR < 51 { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $14 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#           bedtools window -w 10000 -a tmp.bed -b "$ANNOTATION" \
#                > top50.10kb.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 10g,10 -t$'\t' \
#                | awk -F$'\t' 'NR < 101 { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#           bedtools window -w 10000 -a tmp.bed -b "$ANNOTATION" \
#                > top100.10kb.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 10g,10 -t$'\t' \
#                | awk -F$'\t' '{ print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#            bedtools window -w 10000 -a tmp.bed -b "$ANNOTATION" \
#                > all.10kb.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 10g,10 -t$'\t' \
#                | awk -F$'\t' 'NR < 51 { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#            bedtools window -w 0 -a tmp.bed -b "$ANNOTATION" \
#                > top50.within.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 10g,10 -t$'\t' \
#                | awk -F$'\t' 'NR < 101 { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#           bedtools window -w 0 -a tmp.bed -b "$ANNOTATION" \
#                > top100.within.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                | sort -k 10g,10 -t$'\t' \
#                | awk -F$'\t' '{ print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                > tmp.bed \
#                || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#           bedtools window -w 0 -a tmp.bed -b "$ANNOTATION" \
#                > all.within.annotated.tsv \
#                || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#            rm tmp.bed
#
#            cp "$thresh" .
#            "./$(basename "$thresh")" "$comb_output" > empirical_thresholds.txt \
#                || { echo "getting empirical thresholds for {$PHENOTYPE} failed"; exit 1; }
#
#            for threshold in 0.95 0.9; do
#                cutoff="$(grep -P "^${threshold}\t" empirical_thresholds.txt | cut -f 2)"
#                cp "$manhattan" .
#                "./$(basename "$manhattan")" "$FAI" "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" manhattan_0.05.png $cutoff \
#                    || { echo "making Manhattan plot for ${PHENOTYPE} failed"; exit 1; }
#
#                tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.assoc.txt" \
#                    | sort -k 10g,10 -t$'\t' \
#                    | awk -F$'\t' '$10 <'$cutoff' { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
#                    > tmp.bed \
#                    || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
#                bedtools window -w 10000 -a tmp.bed -b "$ANNOTATION" \
#                    > fdr_${threshold}.10kb.annotated.tsv \
#                    || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
#                rm tmp.bed
#
#            done

#            cp *.tsv *.py *.r *.png *.txt *.tsv.gz "${OUTDIR}/${PHENOTYPE}/"
            ;;


        annotate-lm)
            comb_output=combined_output.tsv.gz

            ## Annotation for full set of results
            cd "$SCRATCHDIR/lm/${PHENOTYPE}"
            
            cp "$thresh" .
            "./$(basename "$thresh")" "$comb_output" > empirical_thresholds.txt \
                || { echo "getting empirical thresholds for {$PHENOTYPE} failed"; exit 1; }

            for threshold in 0.95 0.9; do
                cutoff="$(grep -P "^${threshold}\t" empirical_thresholds.txt | cut -f 2)"
                cp "$manhattan" .
                "./$(basename "$manhattan")" "$FAI" "${OUTDIR}/lm/${PHENOTYPE}/output/output.real.effects.assoc.txt" manhattan_0.05.png $cutoff \
                    || { echo "making Manhattan plot for ${PHENOTYPE} failed"; exit 1; }

                tail -n +2 "${OUTDIR}/${PHENOTYPE}/output/output.real.effects.assoc.txt" \
                    | sort -k 10g,10 -t$'\t' \
                    | awk -F$'\t' '$10 <'$cutoff' { print $1 "\t" $3-1 "\t" $3 "\t" $2 ":" $10 };' \
                    > tmp.bed \
                    || { echo "making tmp bed file for ${PHENOTYPE} failed"; exit 1; }
                bedtools window -w 10000 -a tmp.bed -b "$ANNOTATION" \
                    > fdr_${threshold}.10kb.annotated.tsv \
                    || { echo "running bedtools on ${PHENOTYPE} failed"; exit 1; }
                rm tmp.bed
            done

            cp *.tsv *.py *.r *.png *.txt *.tsv.gz "${OUTDIR}/lm/${PHENOTYPE}/"
            ;;


    esac

    # See how long the script really took for future reference
    echo "RUN TIME ${PHENOTYPE} = ${SECONDS} seconds or $(($SECONDS/60)) minutes"

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

    # Run each phenotype in a separate job (of a single array job)
    i=0
    case "$1" in

        prep)
            WALLTIME="00:25:00"
            NTHREADS=1
            MEM=64GB
            echo "$1"$'\t'_$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;

        mvnpermute)
            WALLTIME="00:05:00"
            NTHREADS=1
            MEM=20GB
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        gwas-effects)
            WALLTIME="02:30:00"
            NTHREADS=1
            MEM=8GB
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        gcta)
            WALLTIME="02:00:00"
            NTHREADS=1
            MEM=20GB
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        lm)
            WALLTIME="00:40:00"
            NTHREADS=1
            MEM=8GB
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        gwas-permuted)
            NTHREADS=$NTHREADS_RANDOM
            MEM=96GB
            WALLTIME=56:00:00
            RERUN=no
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE"$'\t'"$RERUN" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE"$'\t'"$RERUN" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        lm-permuted)
            # TODO 2021-07-02: Needs to be re-run
            NTHREADS=$NTHREADS_RANDOM
            MEM=96GB
            WALLTIME=24:00:00
            RERUN=yes
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'"$PHENOFILE"$'\t'"$RERUN" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
#            for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#                echo "$1"$'\t'"$pheno"$'\t'"$PLANT_ONLY_PHENOFILE"$'\t'"$RERUN" > "${ARRAYJOBDIR}/${i}"
#                i=$(($i+1))
#            done
            ;;

        annotate*)
            WALLTIME="04:00:00"
            NTHREADS=1
            MEM=16GB
            for pheno in "${TRAITS[@]}"; do
                echo "$1"$'\t'"$pheno"$'\t'_ > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
           done
#           for pheno in "${PLANT_ONLY_TRAITS[@]}"; do
#               echo "$1"$'\t'"$pheno"$'\t'_ > "${ARRAYJOBDIR}/${i}"
#               i=$(($i+1))
#           done
           ;;
   esac

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi
    array="$(ls "$ARRAYJOBDIR" | tr '\n' ',' | sed 's/,$//g')"

    cd "$WORKDIR"
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --nodes=1 --ntasks=$NTHREADS \
        --time=$WALLTIME --mem=$MEM --array="$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
