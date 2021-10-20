#!/bin/bash
#
# Combine read depth variants with SNPs and small indels called by
# FreeBayes using Mt5.0 as the reference genome. M. truncatula genotypes
# only.
#
# Multi-site variants are broken into allelic primitives.
#
# NOTE: Re-run on 2021-05-26 because I thought that VCF concatenation had a bug,
# then it turned out that the bug doesn't affect this code.
#
# SETTINGS
#
QUEUE="amdsmall"
WALLTIME="09:00:00"
MEM=64GB

PROJDIR="${HOME}/project/select_reseq"
PAVS="${PROJDIR}/../mtr_variant_calls/results/variant_calls/pavs/rdv/2020-11-14/pavs.vcf"
SNPS="${PROJDIR}/../mtr_variant_calls/results/impute/beagle/2020-08-17/imputed.norm.bcf"
STRAINS="${PROJDIR}/data/strain/lists/spanfran2_medicago_truncatula.txt"

OUTDIR="${PROJDIR}/results/combine_variants/spanfran2_medicago_2020-11-16_freebayes_rdvs"
LOGDIR="${OUTDIR}/log"
WORKDIR="${OUTDIR}/working"
SCRIPTDIR="${OUTDIR}/script_copies"

filtertsv="${PROJDIR}/script/bin/filter_tsv_genotypes.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.10.2
    module load htslib/1.9
    module load vcflib/1.0.1
    module load vcftools/0.1.15

    set -euo pipefail

    cd "$OUTDIR"

    cp "$STRAINS" keep.txt

    ## Filter the SNPs:
    ##    keep only individuals in the GWA panel
    ##    remove invariant sites (--min-ac 1:nonmajor)
    ##    already normalized and reduced to primitives
    bcftools view -Ou -S keep.txt "$SNPS" \
        | bcftools view -Ou --min-ac 1:nonmajor \
        > snps.bcf \
        || { echo "processing FreeBayes variants failed"; exit 1; }

    vcfsort "$PAVS" > input_pavs.vcf || { echo "sorting PAVs failed"; exit 1; }
    bgzip -f -i input_pavs.vcf && tabix -p vcf input_pavs.vcf.gz \
        || { echo "tabixing pavs failed"; exit 1; }
    bcftools view -Ov -S keep.txt input_pavs.vcf.gz \
        | vcfallelicprimitives \
        | bcftools norm -Ou -m "-any" \
        | bcftools annotate -Ou --set-id 'rdv-%ID-%CHROM-%POS-%FIRST_ALT' \
        | bcftools view -Ou --min-ac 1:nonmajor \
        > pavs.bcf \
        || { echo "processing read depth variants failed"; exit 1; }
    

    ## Merge
    ## NOTE: This depends on having the same samples in the same order
    ## This was checked manually.
    bcftools view snps.bcf > variants.tmp.vcf \
        || { echo "adding snps failed"; exit 1; }
    bcftools view pavs.bcf \
        | vcfkeepinfo - "." \
        | grep -v "^#" >> variants.tmp.vcf \
        || { echo "adding PAVs failed"; exit 1; }
    bcftools sort variants.tmp.vcf > variants.vcf \
        || { echo "sorting variants failed"; exit 1; }

    ## Compress, index, and tabix
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

    rm input_pavs.* pavs.*cf* snps.bcf *tmp*

    vcftools --gzvcf variants.vcf.gz  --freq2 --out variants \
        || { echo "calculating MAF failed"; exit 1; }

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
        --partition=$QUEUE --nodes=1 --ntasks-per-node=1 \
        --time=$WALLTIME "$sfile"
    echo "$WORKDIR"

fi
