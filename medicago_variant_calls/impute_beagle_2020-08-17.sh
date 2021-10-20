#!/bin/bash
#
# Use BEAGLE 4 to impute variants for M. truncatula accessions.
#
# SETTINGS
#
QUEUE=mangi
WALLTIME=60:00:00
MEM=128GB
NTHREADS=24
#
NITERS=5            # Number of BEAGLE phasing iterations after burnin (default=5); 0 means use old algorithm (I think just burnin)
ERR=0.0001          # Allele error rate (default=0.0001)
NE=1000000          # Effective population size (default=1000000)
MODELSCALE=0.8      # Increasing this reduces runtime & phasing accuracy (default=0.8)
MAXLR=200           # Maximum LR to most likely genotype (default=5000); I think reducing this might improve the speed

PROJDIR="${HOME}/project/mtr_variant_calls"
INFILE="${PROJDIR}/results/variant_calls/snps/freebayes/2020-07-03/variants.qual30.primitives.bcf"
SAMPLES="${PROJDIR}/notes/table/for_imputation_2020-08-17.txt"
FAI="${PROJDIR}/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.fasta.fai"

OUTDIR="${PROJDIR}/results/impute/beagle/2020-08-17"
WORKDIR="${OUTDIR}/working"
ARRAYDIR="${OUTDIR}/arrayjobdata"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"

if [[ "$PBS_ENVIRONMENT" == PBS_BATCH ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load htslib/1.9
    module load bcftools/1.10.2
    module load beagle/20180127

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    REGION="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"

    cd "$OUTDIR"

    case "$TASK" in

        prep)

            ## Convert the bcf input file to vcf.gz
            ## The sample filtering has to take place before the missing data filtering,
            ## then +fixploidy has to be run. If +fixploidy is not run, beagle
            ## won't understand some of the missing genotypes.
            bcftools view -Ou -S "$SAMPLES" "$INFILE" \
                | bcftools view -Ov --include 'N_PASS(GT!="mis")>0' \
                | bcftools +fixploidy -Oz --threads $NTHREADS \
                > tmp.vcf.gz \
                || { echo "converting to vcf.gz format failed"; exit 1; }
            tabix -p vcf tmp.vcf.gz \
                || { echo "tabixing failed"; exit 1; }

            ;;

        beagle)

            ## Run BEAGLE
            java -jar -Xmx"$(echo $MEM | sed 's/GB/g/g')" $BEAGLE_JAR \
                gt=tmp.vcf.gz \
                chrom="$REGION" \
                nthreads=$NTHREADS \
                out="imputed.${REGION}" \
                impute=true \
                gprobs=true \
                niterations=$NITERS \
                ne=$NE \
                err=$ERR \
                modelscale=$MODELSCALE \
                maxlr=$MAXLR \
                || { echo "imputation failed"; exit 1; }

            ;;

        combine)
            ## Combine imputed chromosomes into a single file
            ## Also produce a file with multi-allelic variants split into
            ## multiple records. And filter out sites that are not polymorphic.
            ls imputed.*.vcf.gz | grep -v tmp > to_combine.txt
            while read -r vcf; do
                tabix -f -p vcf "$vcf" \
                    || { echo "tabixxing ${vcf} failed"; exit 1; }
            done < to_combine.txt

            bcftools concat \
                -Ou \
                --threads $NTHREADS \
                -f to_combine.txt \
                | bcftools view \
                -Ob \
                -o imputed.bcf \
                --min-ac 1:nonmajor \
                || { echo "combining imputed variants and dropping non-polymorphic sites failed"; exit 1; }

            bcftools norm -m "-any" -Ou imputed.bcf --threads $NTHREADS \
                | bcftools annotate -Ob --threads $NTHREADS \
                --set-id 'freebayes-mtr-variant-%CHROM-%POS-%FIRST_ALT' \
                > imputed.norm.bcf \
            ;;

            ## Remove intermediate files
            #rm tmp.* *.vcf.gz

    esac

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

    echo "RUN TIME = ${SECONDS} seconds ($(($SECONDS/60)) minutes) for ${REGION}"

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
            WALLTIME=12:00:00
            MEM=20GB
            NTHREADS=$NTHREADS
            echo "$1" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
            ;;
        beagle)
            for chrom in $(cut -f 1 "$FAI"); do
                echo "$1"$'\t'"$chrom" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            done
            ;;

        combine)
            WALLTIME=01:00:00
            MEM=64GB
            NTHREADS=$NTHREADS
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
    qsub2 -A "tiffinp" -W group_list="tiffinp" -q "$QUEUE" \
        -l "walltime=${WALLTIME},nodes=1:ppn=${NTHREADS},mem=${MEM}" \
        -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
