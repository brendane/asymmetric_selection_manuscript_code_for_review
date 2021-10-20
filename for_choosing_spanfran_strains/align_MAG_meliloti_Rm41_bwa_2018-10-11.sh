#!/bin/bash
#
# Use BWA mem to align all the meliloti reads to USDA1106.
# Default BWA settings and uses only paired output from
# the read-cleaning step. Produces cram output.
#
# Also adds read groups. Does not mark duplicates because PCR duplicates
# should have already been removed.
#
# SETTINGS
#
QUEUE="small"
WALLTIME="00:30:00"
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

SPECIES="meliloti"
REFSTRAIN="Rm41_bielefeld"
PROJDIR="/home/tiffinp/epste051/project/jgi_sequencing"
INDIR="${PROJDIR}/results/clean_reads/all/trimgalore_htstream/2018-09-25_dedup_phix_trim/${SPECIES}"
OUTDIR="${PROJDIR}/results/align/MAG_${SPECIES}/${REFSTRAIN}/bwa/2018-10-11"
REFERENCE="${PROJDIR}/data/reference/genome/${REFSTRAIN}/genome.fasta"
WORKDIR="${OUTDIR}/working"
LOGDIR="${OUTDIR}/log"
ARRAYDIR="${OUTDIR}/arrayjobdata"
SCRIPTDIR="${OUTDIR}/script_copies"
SCRATCHDIR="/scratch.global/${USER}/align_jgi_MAG_meliloti_Rm41_bwa_2018-10-11"
depth2cov="${PROJDIR}/script/bin/depth2cov.py"

if [[ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]]; then

    newgrp "tiffinp"
    source "${HOME}/.bashrc"
    source "${HOME}/bin/init-modules.sh"

    module load bwa/0.7.17_gcc-7.2.0_haswell
    module load samtools/1.9

    set -euo pipefail

    STRAIN=$(cut -f 1 -d " " "${ARRAYDIR}/${PBS_ARRAYID}")

    SECONDS=0

    mkdir -p "${OUTDIR}/${STRAIN}"
    cd "$SCRATCHDIR"
    mkdir -p "$STRAIN"
    cd "$STRAIN"

    # Align the paired end reads. Also adds a read group tag in case
    # it will be useful later.
    r1="${INDIR}/${STRAIN}/paired.r1.fastq.gz"
    r2="${INDIR}/${STRAIN}/paired.r2.fastq.gz"
    bwa mem -R "@RG\\tID:${STRAIN}\\tSM:${STRAIN}" \
        -t "$THREADS" -k "$MINSEEDLEN" -w "$BANDWIDTH" -d "$ZDROP" \
        -r "$RESEED" -c "$DISCARDCOUNT" -A "$MATCHSCORE" -B "$MMPENALTY" \
        -O "$GAPOPEN" -E "$GAPXPEN" -L "$CLIPPEN" -U "$UNPPEN" -T "$MAPQ" \
        "$REFERENCE" "$r1" "$r2" | \
        samtools view -b | \
        samtools view -C -F 4 -T "$REFERENCE" | \
        samtools sort -O "CRAM" > "alignment.cram" \
        || { echo "aligning ${STRAIN} failed"; exit 1; }

    # Index the cram file
    samtools index "alignment.cram" \
        || { echo "indexing ${STRAIN} failed"; exit 1; }

    # Get coverage per bp and per replicon
    cp $depth2cov .
    samtools depth -aa "alignment.cram" > "depth.all_sites.tsv" \
        || { echo "getting overall depth for ${STRAIN} failed"; exit 1; }
    cat "depth.all_sites.tsv" | ./$(basename $depth2cov) \
        > "depth.tsv" \
        || { echo "getting mean depth for ${STRAIN} failed"; exit 1; }
    rm -f "depth.all_sites.tsv.gz"
    gzip "depth.all_sites.tsv" \
        || { echo "gzipping depth file for ${STRAIN} failed"; exit 1; }

    # Various coverage and read length statistics
    samtools stats "alignment.cram" > "stats.txt" \
        || { echo "stats failed on ${STRAIN}"; exit 1; }

    # Check for temporary files -- if present, they indicate a failed
    # sorting step that wasn't caught
    ls . | grep "^samtools" && \
        { echo "found samtools temporary file in ${STRAIN} output -- failed"; exit 1; }
   
    # Check for truncated file or invalid format -- just in case of
    # uncaught errors
    samtools flagstat "alignment.cram" | \
        grep -i "invalid\\|truncat" &&
        { echo "invalid cram file for ${STRAIN} failed"; exit 1; }

    cp "alignment.cram" "alignment.cram.crai" "depth.all_sites.tsv.gz" \
        "depth.tsv" "stats.txt" "${OUTDIR}/${STRAIN}/" \
        || { echo "copying files from scratch failed"; exit 1; }

    echo "RUN TIME $SECONDS ($(($SECONDS/60)) minutes) for ${STRAIN}"

    rm "${ARRAYDIR}/${PBS_ARRAYID}"

else

    mkdir -p "$OUTDIR"
    mkdir -p "$LOGDIR"
    mkdir -p "$SCRIPTDIR"
    mkdir -p "$WORKDIR"
    mkdir -p "$ARRAYDIR"
    mkdir -p "$SCRATCHDIR"

    sfile="${SCRIPTDIR}/$(date '+%Y%m%d-%H%M')-$(basename $0)"
    cp "$0" "$sfile"

    i=0
    for strain in "${INDIR}/"*; do 
        if [[ "$strain" =~ "MAG" && ! "$strain" =~ "work" && ! "$strain" =~ "script" && ! "$strain" =~ "log" && ! "$strain" =~ "array" ]]; then
            echo "$(basename "$strain")" > "${ARRAYDIR}/${i}"
            i=$(($i+1))
        fi
    done

    if [[ "$i" == 1 ]]; then
        array="0"
    else
        array="0-$(($i-1))"
    fi

    cd "$WORKDIR"
    qsub2 -A "tiffinp" -W group_list="tiffinp" \
        -q "$QUEUE" -l walltime=$WALLTIME -t "$array" "$sfile"
    echo "$WORKDIR"
    echo "Submitted ${i} jobs"

fi
