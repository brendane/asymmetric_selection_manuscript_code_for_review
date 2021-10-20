#!/bin/bash
#
# Group variants by R^2 at after filtering on MAF and missing
# data. PAVs called from SibeliaZ-based whole genome alignment,
# with separate alignments for medicae and meliloti, and SNPs from FreeBayes.
#
# SETTINGS
#
WALLTIME=01:00:00
NTHREADS=1
MEM=8GB
QUEUE=amdsmall
#
MIN_MAF=0.05
R2S=(0.20 0.80 0.90 0.95 0.99)
MAX_MISS=0.20
# SP1-L = SpanFran1 (summer 2019) 89 strain meliloti (MAG194-mel dropped from analysis to make 88 strains),
# SP1-D = SpanFran1 89 strain medicae (MAG194-med dropped from analysis to make 88 strains),
# C101 = original 101 strain meliloti (PNAS + legacy),
# C68 = 68 strain (NPD & winter 2018) meliloti,
# SP2 = SpanFran2 (big GWAS) 88 strain meliloti

PROJDIR="${HOME}/project/select_reseq"
VARIANTS_MEL="${PROJDIR}/results/combine_variants/meliloti_2021-04-06_denovo/variants.bcf"
VARIANTS_MED="${PROJDIR}/results/combine_variants/medicae_2021-04-06_denovo/variants.bcf"
OLD_PYTHONPATH="$HOME/usr/lib/python2.7/site-packages.old"

declare -A STRAINS
STRAINS[SP2]="${PROJDIR}/data/strain/lists/spanfran2_meliloti.txt"
STRAINS[SP-C86]="${PROJDIR}/data/strain/lists/spanfran_c86.txt"
STRAINS[SP1-L]="${PROJDIR}/data/strain/lists/spanfran1_meliloti.txt"
STRAINS[C68]="${PROJDIR}/data/strain/lists/68_strains.txt"
STRAINS[C101]="${PROJDIR}/data/strain/lists/101Strains.txt"
STRAINS[SP1-D]="${PROJDIR}/data/strain/lists/spanfran1_medicae.txt"

declare -A VARIANTS
VARIANTS[SP1-L]="$VARIANTS_MEL"
VARIANTS[C68]="$VARIANTS_MEL"
VARIANTS[C101]="$VARIANTS_MEL"
VARIANTS[SP2]="$VARIANTS_MEL"
VARIANTS[SP-C86]="$VARIANTS_MEL"
VARIANTS[SP1-D]="$VARIANTS_MED"

OUTDIR="${PROJDIR}/results/ld/pangenome/r2_groups/2021-04-06_denovo"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"

vcf2tsv="${PROJDIR}/script/bin/vcf2tsv.py"
grouper="${PROJDIR}/script/bin/r2_groups_xtra_sort"
chooser="${PROJDIR}/script/bin/ld_grouping_choose_one.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.10.2

    set -euo pipefail

    SECONDS=0

    COMMUNITY="$(cut -f 1 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"
    R2="$(cut -f 2 "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}")"

    mkdir -p "${OUTDIR}/${COMMUNITY}/${R2}" && cd "${OUTDIR}/${COMMUNITY}/${R2}"

    case "$COMMUNITY" in
        SP*)
            cat "${STRAINS[$COMMUNITY]}" \
                | grep -v "MAG194" > strains.txt \
                || { echo "making list of strain names to match variants file failed"; exit 1; }
            ;;
        C101)
            sed 's/^/USDA/g' "${STRAINS[$COMMUNITY]}" \
                | sed 's/USDA3011/3011/g' > strains.txt \
                || { echo "making list of strain names to match variants file failed"; exit 1; }
            ;;
        C68)
            awk '/^[A-Z]/ {print $1}; /^[0-9]/ {print "USDA" $1}' \
                "${STRAINS[$COMMUNITY]}" \
                | sed 's/USDA3085/3085/g' \
                > strains.txt \
                || { echo "making list of strain names to match variants file failed"; exit 1; }
            ;;
    esac

    ## Filter; make sure info field is removed so that bcftools
    ## will calculate MAF and missingness using the strains in the
    ## file.
    bcftools view -Ou -S strains.txt "${VARIANTS[$COMMUNITY]}" \
        | bcftools view -Ou --min-af "${MIN_MAF}:minor" \
        | bcftools view -Ov -e 'F_PASS(GT=="mis")>'"${MAX_MISS}" \
        > filtered_variants.vcf \
        || { echo "filtering variants failed"; exit 1; }

    cp "$vcf2tsv" .
    "./$(basename "$vcf2tsv")" --rs --output filtered_variants.tsv \
        filtered_variants.vcf \
        || { echo "converting variants to tsv failed"; exit 1; }

    n_strains=$(wc -l strains.txt | cut -f 1 -d" ")
    n_strains_tsv=$(awk 'NR==1 { print NF-3 };' filtered_variants.tsv)
    if [[ $n_strains != $n_strains_tsv ]]; then
        echo "not all strains were found, failed"
        exit 1
    fi

    # Do the grouping
    cp "$grouper" .
    ./$(basename "$grouper") $R2 1 filtered_variants.tsv \
        > tmp.groups.tsv \
        || { echo "grouping at ${R2} in ${COMMUNITY} failed"; exit 1; }

    # Gzip the variants file for later use
    gzip -f filtered_variants.tsv \
        || { echo "gzipping variants for ${COMMUNITY} failed"; exit 1; }

    # Sort by group and by seed status
    head -n 1 tmp.groups.tsv > output.tsv
    tail -n +2 tmp.groups.tsv | sort -t$'\t' -k 6n,6 -k 7rn,7 >> output.tsv \
        || { echo "sorting at ${R2} in ${COMMUNITY} failed"; exit 1; }

    # Choose one variant to represent each group; given the way
    # the file is sorted, this will be the variant with the greatest
    # number of individuals genotyped, and if there are ties, the
    # greatest minor allele frequency.
    cp "$chooser" .
    PYTHONPATH="$OLD_PYTHONPATH" ./$(basename "$chooser") \
        --output one_variant.tsv output.tsv \
        || { echo "choosing one variant for ${R2} in ${COMMUNITY} failed"; exit 1; }

    # Eliminate temporary files
#    rm tmp.* filtered_variants.vcf

    echo "RUN TIME = ${SECONDS} ($(($SECONDS/60)) minutes)"

    rm "${ARRAYDIR}/${SLURM_ARRAY_TASK_ID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y-%m-%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    for community in "${!STRAINS[@]}"; do
        for R2 in "${R2S[@]}"; do
            echo "$community"$'\t'"$R2" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
        done
    done

    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi

    cd $WORKDIR
    sbatch --account=tiffinp --export=PBS_ENVIRONMENT=PBS_BATCH \
        --partition=$QUEUE --mem=$MEM --nodes=1 \
        --ntasks-per-node=$NTHREADS \
        --time=$WALLTIME --array="$array" "$sfile"
    echo $WORKDIR
    echo "Submitted ${i} jobs"

fi
