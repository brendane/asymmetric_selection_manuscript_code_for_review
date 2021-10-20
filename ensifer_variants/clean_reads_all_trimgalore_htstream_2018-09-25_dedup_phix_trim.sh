#!/bin/bash
#
# Clean the reads by removing PCR duplicates, adaptors, and PhiX
# sequences, and trimming low quality bases.
#
# Data:
#   All medicae and meliloti raw reads from JGI as of 25 Sept 2018 and
#   all previous medicae and meliloti sequences (2010 and 2016
#   sequencing) downloaded from NCBI.
#
# Steps:
#   1. PCR duplicate removal with HT-Stream
#   2. Adaptor trimming with HT-Stream
#   3. PhiX removal with HT-Stream
#   4. Quality trimming and more adaptor removal with TrimGalore!
#
# Runs on MSI.
#

QUEUE="small"
MEM="16GB"
WALLTIME="05:00:00"

# HT-Stream settings
DEDUP_LOG_FREQ=0  # Turns off logging of duplicates
DEDUP_MINLEN=20   # Length of unique ID (default=10)
DEDUP_AVG_Q=35    # Minimum average quality score to have the read written automatically if a duplicate has not been found yet (default=30)
DEDUP_IAVG_Q=10   # Min. avg. qual. score to consider a read informative (default=5)
DEDUP_START=10    # Starting point in read of unique ID (default=10)
ATRIM_FIXBASES="--no-fixbases" # Do not change overlapping bases to consensus
ATRIM_MMEDENS="0.1" # Maximum proportion (percentage?) of mismatches in overlapped section (default=0.25)
ATRIM_ORPHANS="--no-orphans" # Do not write out orphaned SE reads

# TrimGalore! settings
# Same settings for all types of reads, regardless of original length
QUALTYPE="sanger"
QUALTHRESH=30
LENTHRESH=80
SLENTHRESH=81
ADAPTSTRINGENCY=3
ERRORRATE=0.1

PROJDIR="${HOME}/project/jgi_sequencing"
READSDIR="${PROJDIR}/data/reads"
INDIRS=("${READSDIR}/ncbi/medicae" "${READSDIR}/ncbi/meliloti"
        "${READSDIR}/jgi/heath_medicae" "${READSDIR}/jgi/heath_meliloti")
OUTDIR="${PROJDIR}/results/clean_reads/all/trimgalore_htstream/2018-09-25_dedup_phix_trim"

WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/jgi-read-processing-2018-09-25"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp tiffinp
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load htslib/1.6
    module load htstream/1.0.0
    module load trimgalore/0.5.0

    set -euo pipefail

    SECONDS=0

    R1_raw="$(cut -f 1 "${ARRAYDIR}/${PBS_ARRAYID}")"
    SPP="$(cut -f 2 "${ARRAYDIR}/${PBS_ARRAYID}")"
    R2_raw="${R1_raw/_R1.fastq/_R2.fastq}"
    STRAIN="$(basename "${R1_raw/_R1.fastq.gz/}")"

    cd "$SCRATCHDIR"
    mkdir -p "$SPP" && cd "$SPP" && mkdir -p "$STRAIN" && cd "$STRAIN"

    ## Added after processing all samples except for E. medicae M161 (an NCBI sample)
    cp "$R1_raw" "$R2_raw" .
    gunzip -f "$(basename "$R1_raw")"
    gunzip -f "$(basename "$R2_raw")"
    R1="$(basename "${R1_raw/.gz}")"
    R2="$(basename "${R2_raw/.gz}")"

    # 1. Remove PCR duplicates
    # 2. Adaptor trimming by looking at read overlaps
    # 3. PhiX contaminant removal
    hts_SuperDeduper -1 "$R1" -2 "$R2" --tab-output --to-stdout \
        --log_freq "$DEDUP_LOG_FREQ" --length "$DEDUP_MINLEN" \
        --avg-qual-score "$DEDUP_AVG_Q" \
        --inform-avg-qual-score "$DEDUP_IAVG_Q" --start "$DEDUP_START" \
        --stats-file "${STRAIN}.htseq.stats.txt" | \
        hts_AdapterTrimmer --from-stdin \
            "$ATRIM_FIXBASES" "$ATRIM_ORPHANS" \
            --max-mismatch-errorDensity "$ATRIM_MMEDENS" \
            --tab-output --to-stdout \
            --stats-file "${STRAIN}.htseq.stats.txt" --append-stats-file | \
        hts_SeqScreener --from-stdin \
            --fastq-output -p "${STRAIN}.htseq" --force \
            --stats-file "${STRAIN}.htseq.stats.txt" --append-stats-file \
        || { echo "HT-Stream pipeline failed on ${STRAIN}"; exit 1; }

    # 4. Further adaptor removal and quality trimming with TrimGalore!
    htst1="${STRAIN}.htseq_R1.fastq"
    htst2="${STRAIN}.htseq_R2.fastq"
    o1="${htst1/.fastq/_val_1.fq.gz}"
    o2="${htst2/.fastq/_val_2.fq.gz}"
    u1="${htst1/.fastq/_unpaired_1.fq.gz}"
    u2="${htst2/.fastq/_unpaired_2.fq.gz}"
    trim_galore --quality "$QUALTHRESH" --stringency "$ADAPTSTRINGENCY" \
        -e "$ERRORRATE" --gzip --length "$LENTHRESH" --paired \
        --retain_unpaired --length_1 "$SLENTHRESH" --length_2 "$SLENTHRESH" \
        "$htst1" "$htst2" | tee "log.txt" \
        || { echo "running trim galore on ${STRAIN} failed"; exit 1; }
    echo "cutadapt version:" $(cutadapt --version) >> "log.txt"
    zcat "$u1" "$u2" | bgzip -c > \
        "single.fastq.gz" \
        || { echo "concatenating singles for ${STRAIN} failed"; exit 1; }
    rm "$u1" "$u2"
    mv "$o1" "paired.r1.fastq.gz" \
        || { echo "moving R1 output for ${STRAIN} failed"; exit 1; }
    mv "$o2" "paired.r2.fastq.gz" \
        || { echo "moving R1 output for ${STRAIN} failed"; exit 1; }

    # 5. Count reads
    rm -f "read_counts.txt"
    echo "pairs	$(zcat "paired.r1.fastq.gz" | blwc - | sed 's/pairs.\\+\\t/pairs/g')" > "read_counts.txt" \
        || { echo "counting number of read pairs left in ${STRAIN} failed"; exit 1; }
    echo "singles	$(zcat "single.fastq.gz" | blwc - | sed 's/singles.\\+\\t/singles/g')" >> "read_counts.txt" \
        || { echo "counting number of single read left in ${STRAIN} failed"; exit 1; }

    # 6. Copy files to output directory; leaving out the unpaired reads
    #    file b/c it hasn't been used lately.
    mkdir -p "${OUTDIR}/${SPP}/${STRAIN}"
    cp "read_counts.txt" "paired.r1.fastq.gz" \
        "paired.r2.fastq.gz" "log.txt" "${STRAIN}.htseq.stats.txt" \
        "${OUTDIR}/${SPP}/${STRAIN}" \
        || { echo "copying files for ${STRAIN} failed"; exit 1; }

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

    echo "RUN TIME FOR ${STRAIN} = ${SECONDS} ($(($SECONDS/60)) minutes)"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(basename "$0")-$(date "+%Y-%d-%m-%H%M")"
    cp "$0" "$sfile"

    cd "$WORKDIR"
    i=0
    # Note that this will not work if there are or newlines tabs in the
    # file names, but these files are safe.
    for id in "${INDIRS[@]}"; do
        for r1 in "$id"/*_R1.fastq.gz; do
            # This if statement was added temporaryily:
            if [[ "$r1" =~ "USDA1150" || "$r1" =~ "USDA1169" ]]; then
                species="meliloti"
                if [[ "$r1" =~ "medicae/" ]]; then
                    species="medicae"
                fi
                echo "${r1}"$'\t'"${species}" > "${ARRAYDIR}/${i}"
                i=$(($i+1))
            fi
        done
    done
    if [[ "$i" == 1 ]]; then
        array=0
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub2 -q "$QUEUE" -l "walltime=${WALLTIME},mem=${MEM}" -t "$array" "$sfile"
    echo "Submitted ${i} jobs"
    echo "$WORKDIR"

fi
