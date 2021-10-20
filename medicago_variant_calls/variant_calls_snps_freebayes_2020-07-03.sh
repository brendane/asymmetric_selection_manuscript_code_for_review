#!/bin/bash
#
# Use FreeBayes to call SNPs and small indels on all the Mtr accessions.
#
# SETTINGS
#
QUEUE=mangi
WALLTIME=25:00:00
MEM=20GB
#MEM=60GB # One region needed needed extra memory
NTHREADS=1
NTHREADS_BCFTOOLS=4
#
REGIONSIZE=300000 # Size of regions to split the reference into for parallel operation
MAPQ=30           # Minimum mapping quality required for reads
POP_PRIORS="--no-population-priors" # Do not call variants based on population allele frequencies
N_ALLELES=6       # No more than this many alleles per variant (to save memory)
PLOIDY=2
VARIANT_QUAL_THRESHOLD=30 # Minimum variant quality for final output file

PROJDIR="${HOME}/project/mtr_variant_calls"
INDIR="${PROJDIR}/results/align/mt5/bwa/2020-05-08"
OUTDIR="${PROJDIR}/results/variant_calls/snps/freebayes/2020-07-03"
REFERENCE="${PROJDIR}/data/assemblies/Mtr_genomes/MtrunA17r5.0-ANR/genome.fasta"

WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
SCRIPTDIR="${OUTDIR}/script_copies"
ARRAYJOBDIR="${OUTDIR}/arraydata"
SCRATCHDIR="/scratch.global/${USER}/variant_calls_snps_MTR_freebayes_2020-07-03"
TEMPVCFDIR="${SCRATCHDIR}/tmpvcf"

if [[ "$PBS_ENVIRONMENT" == PBS_BATCH ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bcftools/1.9
    module load freebayes/20200703
    module load htslib/1.9
    module load vcflib/1.0.1

    SECONDS=0

    TASK="$(cut -f 1 "${ARRAYJOBDIR}/${PBS_ARRAYID}")"
    REGION="$(cut -f 2 "${ARRAYJOBDIR}/${PBS_ARRAYID}")"

    cd "$SCRATCHDIR"

    case "$TASK" in

        "prep")
            rm -f "local_bams.txt"
            IFS=$'\n' ls "$INDIR/" \
                | grep -v "work\|array\|log\|script\|\." \
                > "accessions.txt" \
                || { echo "listing strains failed"; exit 1; }

            while read -r geno; do
                for crfile in "${INDIR}/${geno}/"*.cram; do
                    local_name="$(pwd)/${geno}-$(basename "$crfile")"
                    cp "$crfile" "$local_name" \
                        || { echo "copying ${crfile} cram file failed"; exit 1; }
                    cp "${crfile}.crai" "${local_name}.crai" \
                        || { echo "copying ${crfile} crai file failed"; exit 1; }
                    echo "$local_name" >> "local_bams.txt"
                done
            done < "accessions.txt"
            ;;

        "run")
            mkdir -p "$TEMPVCFDIR"
            rm -f "${TEMPVCFDIR}/tmp.${REGION}.vcf"
            freebayes \
                --region "$REGION" \
                -f "$REFERENCE" \
                -p "$PLOIDY" \
                -L "local_bams.txt" \
                --min-mapping-quality "$MAPQ" \
                $POP_PRIORS \
                --use-best-n-alleles "$N_ALLELES" \
                --gvcf \
                > "${TEMPVCFDIR}/tmp.${REGION}.vcf" \
                || { echo "calling variants for ${REGION} failed"; exit 1; }
            ;;

        "combine")
#            ## Concatenate all the tmp files using the same strategy
#            ## employed by freebayes-parallel
#            cat $(ls "$TEMPVCFDIR/"* | sort -t: -k 1b,1 -k 2n,2) \
#                | vcffirstheader \
#                | vcfstreamsort -w 1000 \
#                | vcfuniq \
#                > "merged.gvcf" \
#                || { echo "combining vcfs failed"; exit 1; }
#
#            ## Pick out only putative variant sites from the gvcf by
#            ## setting a low minimum quality threshold.
#            vcffilter -f "QUAL > 1" merged.gvcf > tmp.vcf \
#                || { echo "picking putative variant sites failed"; exit 1; }
#
#            ## Prepare for use with bcftools
#            rm -f tmp.vcf.gz
#            bgzip -i tmp.vcf \
#                || { echo "gzipping tmp file failed"; exit 1; }
#            tabix tmp.vcf.gz \
#                || { echo "tabixing tmp file failed"; exit 1; }
#
#            # Add variant IDs, change name, and compress and index
#            bcftools annotate --set-id 'freebayes-snp-%CHROM-%POS-%FIRST_ALT' tmp.vcf.gz \
#                > variants.vcf \
#                || { echo "normalizing SNPs failed"; exit 1; }
#            rm -f variants.vcf.gz
#            bgzip -i variants.vcf \
#                || { echo "gzipping variants failed"; exit 1; }
#            tabix variants.vcf.gz \
#                || { echo "tabixing variants failed"; exit 1; }
#
#            # Filter variants by quality score, then compress and index
#            vcffilter -f "QUAL > ${VARIANT_QUAL_THRESHOLD}" variants.vcf.gz \
#                > variants.qual${VARIANT_QUAL_THRESHOLD}.vcf \
#                || { echo "filtering at QUAL >= ${VARIANT_QUAL_THRESHOLD} failed"; exit 1; }
#            rm -f variants.qual${VARIANT_QUAL_THRESHOLD}.vcf.gz
#            bgzip -i variants.qual${VARIANT_QUAL_THRESHOLD}.vcf \
#                || { echo "gzipping filtered variants failed"; exit 1; }
#            tabix variants.qual${VARIANT_QUAL_THRESHOLD}.vcf.gz \
#                || { echo "tabixing filtered variants failed"; exit 1; }
#
#            bcftools view --threads $NTHREADS_BCFTOOLS \
#                -O b -o "variants.qual${VARIANT_QUAL_THRESHOLD}.bcf" \
#                "variants.qual${VARIANT_QUAL_THRESHOLD}.vcf.gz" \
#                || { echo "converting to bcf format failed"; exit 1; }
#
#            vcfallelicprimitives "variants.qual${VARIANT_QUAL_THRESHOLD}.vcf.gz" \
#                | bcftools norm -O u -f "$REFERENCE" -d any \
#                | bcftools annotate -O u --set-id 'freebayes-var-%CHROM-%POS-%FIRST_ALT' \
#                | bcftools view --threads $NTHREADS_BCFTOOLS -O b \
#                -o "variants.qual${VARIANT_QUAL_THRESHOLD}.primitives.bcf" \
#                || { echo "making allelic primitives file and converting to bcf failed"; exit 1; }
#            bcftools index "variants.qual${VARIANT_QUAL_THRESHOLD}.primitives.bcf" \
#                || { echo "indexing primitives file failed"; exit 1; }
#
#            # Statistics at different quality filter stringencies
#            vcfstats variants.vcf.gz > stats.all.txt \
#                || { echo "getting stats failed"; exit 1; }
#            for q in 30; do
#                vcffilter -f "QUAL > ${q}" variants.vcf.gz \
#                    | vcfstats \
#                    > "stats.qual${q}.txt" \
#                    || { echo "getting stats for qual > ${q} failed"; exit 1; }
#            done

            # Copy from scratch directory
            # Note that gvcf should go into second tier storage because
            # it will very large
            cp variants.qual${VARIANT_QUAL_THRESHOLD}.bcf "$OUTDIR" \
                || { echo "copying vcf and indexes failed"; exit 1; }
            cp variants.qual${VARIANT_QUAL_THRESHOLD}.primitives.bcf* "$OUTDIR" \
                || { echo "copying vcf and indexes failed"; exit 1; }
            cp stats.* "$OUTDIR" \
                || { echo "copying stats files failed"; exit 1; }

            ## Gzip the merged.gvcf
            rm -f merged.gvcf.gz
            bgzip -i merged.gvcf \
                || { echo "gzipping gvcf failed"; exit 1; }
            tabix merged.gvcf.gz \
                || { echo "tabixing gvcf failed"; exit 1; }

    esac

    rm "${ARRAYJOBDIR}/${PBS_ARRAYID}"

    echo "RUN TIME = ${SECONDS} seconds ($(($SECONDS/60)) minutes)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYJOBDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    case "$1" in
        "prep")
            WALLTIME=01:00:00
            MEM=2GB
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
        "run")
            module load freebayes/20200703
            #rm -rf "$TEMPVCFDIR"
            for region in $(fasta_generate_regions.py "${REFERENCE}.fai" "$REGIONSIZE"); do
                echo "$1"$'\t'"$region" > "${ARRAYJOBDIR}/${i}"
                i=$(($i+1))
            done
            ;;
        "combine")
            NTHREADS=$NTHREADS_BCFTOOLS
            WALLTIME=72:00:00
            MEM=64GB
            echo "$1"$'\t'"_" > "${ARRAYJOBDIR}/${i}"
            i=$(($i+1))
            ;;
    esac

    if [[ "$i" == 1 ]]; then
        array=0
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
