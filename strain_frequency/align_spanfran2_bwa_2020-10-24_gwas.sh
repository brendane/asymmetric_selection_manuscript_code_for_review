#!/bin/bash
#
# Use BWA mem to align the winter 2020 SpanFran2 GWAS experiment
# reads to the reference genome.
# Default BWA settings and uses only paired output from
# TrimGalore!. Produces cram output.
#
# Also adds read groups and marks duplicates.
#
# Starts with reads trimmed with TrimGalore!.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="05:00:00"
#
# bwa settings - all defaults; placed here for record-keeping
#
THREADS=1           # -t option
MINSEEDLEN=19       # -k option
BANDWIDTH=100       # -w option
ZDROP=100           # -d option
RESEED=1.5          # -r option
DISCARDCOUNT=10000  # -c option
MATCHSCORE=1        # -A option
MMPENALTY=4         # -B option
GAPOPEN=6           # -O option
GAPXPEN=1           # -E option
CLIPPEN=5           # -L option
UNPPEN=9            # -U option
MAPQ=30             # -T option (min alignment score)

PROJDIR="${HOME}/project/select_reseq"
DATASET="spanfran2"
INDIR="${PROJDIR}/results/trim_reads/${DATASET}/trim_galore/2020-10-23_gwas"
REFERENCE="${PROJDIR}/data/reference/USDA1106/USDA1106.final.fasta"

OUTDIR="${PROJDIR}/results/align/${DATASET}/bwa/2020-10-24_gwas"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/align_spanfran2_2020-10-24"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bwa/0.7.17_gcc-7.2.0_haswell
    module load samtools/1.9

    set -euo pipefail

    POOL=$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")

    SECONDS=0

    mkdir -p "${OUTDIR}/${POOL}"
    cd "$SCRATCHDIR"
    mkdir -p "$POOL" && cd "$POOL"

    # Align the paired end reads. Also adds a read group tag in case
    # it will be useful later.
    r1="${INDIR}/${POOL}/paired.r1.fastq.gz"
    r2="${INDIR}/${POOL}/paired.r2.fastq.gz"
    bwa mem -R "@RG\\tID:${DATASET}_${POOL}\\tSM:${POOL}" \
        -t "$THREADS" -k "$MINSEEDLEN" -w "$BANDWIDTH" -d "$ZDROP" \
        -r "$RESEED" -c "$DISCARDCOUNT" -A "$MATCHSCORE" -B "$MMPENALTY" \
        -O "$GAPOPEN" -E "$GAPXPEN" -L "$CLIPPEN" -U "$UNPPEN" -T "$MAPQ" \
        "$REFERENCE" "$r1" "$r2" | \
        samtools view -b | \
        samtools sort -n > "initial.bam" \
        || { echo "aligning ${POOL} failed"; exit 1; }

    # Mark duplicates
    samtools fixmate --reference "$REFERENCE" -m "initial.bam" - | \
        samtools sort - | \
        samtools markdup --reference "$REFERENCE" - - | \
        samtools sort --reference "$REFERENCE" -O "CRAM" - > "alignment.cram" \
        || { echo "marking duplications failed"; exit 1; }

    # Index the cram file
    samtools index "alignment.cram" \
        || { echo "indexing ${POOL} failed"; exit 1; }

    # Various coverage and read length statistics
    samtools stats "alignment.cram" > "stats.txt" \
        || { echo "stats failed on ${POOL}"; exit 1; }

    # Check for temporary files -- if present, they indicate a failed
    # sorting step that wasn't caught
    # NOTE: should add -a to ls command
    ls . | grep "^samtools" && \
        { echo "found samtools temporary file in ${POOL} output -- failed"; exit 1; }
   
    # Check for truncated file or invalid format -- just in case of
    # uncaught errors
    samtools flagstat "alignment.cram" | \
        grep -i "invalid\\|truncat" &&
        { echo "invalid cram file for ${POOL} failed"; exit 1; }

    cp "alignment.cram" "alignment.cram.crai" \
        "stats.txt" "${OUTDIR}/${POOL}/" \
        || { echo "copying files from scratch failed"; exit 1; }

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes) for ${POOL}"

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp $0 $sfile

    i=0
    for pool in $(ls "$INDIR" | grep -v "work\|script\|log\|array"); do
        echo $pool > "${ARRAYDIR}/${i}"
        i=$(($i+1))
    done

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q $QUEUE -l walltime=$WALLTIME -t "$array" $sfile
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
